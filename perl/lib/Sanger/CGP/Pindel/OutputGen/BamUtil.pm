package Sanger::CGP::Pindel::OutputGen::BamUtil;

########## LICENCE ##########
# Copyright (c) 2014-2018 Genome Research Ltd.
#
# Author: CASM/Cancer IT <cgphelp@sanger.ac.uk>
#
# This file is part of cgpPindel.
#
# cgpPindel is free software: you can redistribute it and/or modify it under
# the terms of the GNU Affero General Public License as published by the Free
# Software Foundation; either version 3 of the License, or (at your option) any
# later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
# details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
########## LICENCE ##########


use Sanger::CGP::Pindel;
use Sanger::CGP::Pindel::Implement;
use Sanger::CGP::Vcf::Contig;
use Sanger::CGP::Vcf::Sample;

use strict;
use autodie qw(:all);
use Carp;
use Data::Dumper;
use File::Which qw(which);
use Bio::DB::HTS;
use PerlIO::gzip;

use Const::Fast qw(const);

const my $PG_TEMPLATE => "\@PG\tID:%s\tPN:%s\tCL:%s\tPP:%s\tDS:%s\tVN:%s";
const my $BAMSORT => ' inputformat=%s outputformat=%s reference=%s | tee %s | md5sum -b - > %s';

sub sam_to_sorted_bam {
  my ($path_prefix, $base_dir, $sam_files, $as_cram, $as_csi, $ref) = @_;
  $sam_files = Sanger::CGP::Pindel::Implement::fragmented_files($base_dir, $sam_files, '@', 'FINAL_MERGED.sam.gz');
  update_header_when_no_reads($base_dir, $sam_files);
  my $aln_fmt = 'bam';
  my $idx_fmt = 'bai';
  my $idx_switch = q{-b};
  if($as_cram) {
    $aln_fmt = 'cram';
    $idx_fmt = 'crai';
    $idx_switch = q{};
  }
  elsif($as_csi) {
    $idx_fmt = 'csi';
    $idx_switch = q{-c};
  }
  my $aln_out = $path_prefix.'.'.$aln_fmt;
  my $aln_idx = $aln_out.'.'.$idx_fmt;
  my $aln_md5 = $aln_out.'.md5';
  my $command = q{};
  if($sam_files->[0] =~ m/FINAL_MERGED[.]sam.gz$/) {
    $command .= sprintf ' zcat %s | ', $sam_files->[0];
  }
  else {
    my $header_from = $sam_files->[0];
    $command .= "cd $base_dir; ";
    $command .= sprintf q{ ( zcat %s | samtools view -H - ; zcat %s | grep -v '^@' || true ) | }, $header_from, join q{ }, @{$sam_files};
  }

  $command .= which('bamsort');
  $command .= sprintf $BAMSORT, 'sam', $aln_fmt, $ref, $aln_out, $aln_md5;

  $command = sprintf q{bash -c "set -o pipefail; %s"}, $command;
  warn "running: $command";
  system($command);

  $command = sprintf 'samtools index %s %s %s', $idx_switch, $aln_out, $aln_idx;
  warn "running: $command";
  system($command);

  return 1;
}

sub update_header_when_no_reads {
  my ($base_dir, $sam_files) = @_;
  for my $sam_file (@{$sam_files}) {
    my $sam;
    if($sam_file =~ m/FINAL_MERGED[.]sam.gz$/) {
      $sam = $sam_file
    }
    else {
      $sam = $base_dir.$sam_file;
    }
    my @header;
    my $has_records = 0;
    open my $IN, '<:gzip', $sam || die $!;
    while(<$IN>) {
      if($_ =~ m/^\@/) {
        chomp $_;
        push @header, $_;
      }
      else {
        $has_records++;
      }
      last if($has_records);
    }
    close $IN;
    unless($has_records) {
      if($header[0] =~ s/SO:unknown/SO:coordinate/) {
        warn "Updating sort order as no records\n";
        open my $OUT, '>:gzip', $sam || die $!;
        print $OUT join("\n",@header),"\n";
        close $OUT;
      }
    }
  }
}

sub pindel_header {
  my ($bam, $samp_type, $parent_pg, $cmd, $man_species, $man_assembly) = @_;

  my $hfile = Bio::DB::HTSfile->open($bam);
  my @h_lines = split /\n/, $hfile->header_read->text();
  my $last_pg;
  for my $idx(0..(scalar @h_lines)-1) {
    my $h_line = $h_lines[$idx];
    if($h_line =~ m/^\@SQ/) {
      $h_lines[$idx] .= "\tSP:$man_species" if(defined $man_species && $h_line !~ m/\tSP:/);
      $h_lines[$idx] .= "\tAS:$man_assembly" if(defined $man_assembly && $h_line !~ m/\tAS:/);
    }
    $last_pg = $idx if($h_line =~ m/^\@PG/);
  }

  if(defined $last_pg) {
    my ($last_id) = $h_lines[$last_pg] =~ m/\tID:([^\t]+)/;
    my $preceeding_id;
    if($parent_pg) {
      ($preceeding_id) = $parent_pg =~ m|\\tID:([^\\]+)\\t|;
      $parent_pg =~ s/\\t/\t/g;
      $parent_pg =~ s/\tPP:\./\tPP:$last_id/;
      push @h_lines, $parent_pg;
      ($last_id) = $h_lines[-1] =~ m/\tID:[^\t]+/;
    }
    push @h_lines, pg_from_caller('pindel_to_combined_vcf', 'Converts text output of pindel to VCF and BAM alinments', $VERSION, $cmd, $preceeding_id);
  }
  return join "\n", @h_lines;
}

sub pg_from_caller {
  my ($id, $desc, $version, $cmd, $pp_id) = @_;
  $cmd =~ s/\t/\\t/g;

  return sprintf $PG_TEMPLATE,  $id,
                                $0, # program name
                                $cmd,
                                $pp_id || q{.},
                                $desc,
                                $version;
}

sub parse_contigs{
  my($header_txt, $man_species, $man_assembly) = @_;

  my $contigs = {};

  foreach my $line (split(/\n/,$header_txt)){
    if($line =~ /^\@SQ/){
      my ($name) = $line =~ /SN:([^\t]+)/;
      my ($length) = $line =~ /LN:([^\t]+)/;
      my ($assembly) = $line =~ /AS:([^\t]+)/;
      my ($species) = $line =~ /SP:([^\t]+)/;

      if(defined $man_species) {
        if(defined $species) {
          warn "WARN: Manually entered species ($man_species) doesn't match BAM file header ($species). Defaulting to BAM header species.\n"
            if($man_species ne $species);
        }
        else { $species = $man_species; }
      }
      die "ERROR: Species must be defined, check options/BAM header\n" unless(defined $species);

      if(defined $man_assembly) {
        if(defined $assembly) {
          warn "WARN: Manually entered assembly ($man_assembly) doesn't match BAM file header ($assembly). Defaulting to BAM header assembly.\n"
            if($man_assembly ne $assembly);
        }
        else { $assembly = $man_assembly; }
      }
      die "ERROR: Assembly must be defined, check options/BAM header\n" unless(defined $assembly);

      my $contig = new Sanger::CGP::Vcf::Contig(
        -name => $name,
        -length => $length,
        -assembly => $assembly,
        -species => $species
      );

      if(exists $contigs->{$name}){
      	croak "ERROR: Trying to merge contigs with conflicting data:\n".Dumper($contigs->{$name})."\n".Dumper($contig)
          unless _compare($contig, $contigs->{$name});
      } else {
      	$contigs->{$name} = $contig;
      }
    }
  }
  return $contigs;
}

sub parse_samples{
  my($header_txt,$study,$protocol,$accession,$accession_source,$description) = @_;

  my $samples = {};

  foreach my $line (split(/\n/,$header_txt)){
    if($line =~ /^\@RG/){

      my ($name) = $line =~ /SM:([^\t]+)/;
      my ($platform) = $line =~ /PL:([^\t]+)/;

      $samples->{$name} = new Sanger::CGP::Vcf::Sample(
        -name => $name,
        -study => $study,
        -platform => $platform,
        -seq_protocol => $protocol,
        -accession => $accession,
        -accession_source => $accession_source,
        -description => $description
      ) unless exists $samples->{$name};
    }
  }
  return $samples;
}


sub _compare{
	my($contig,$contig1) = @_;
	croak 'Incorrect input' unless ref $contig eq 'Sanger::CGP::Vcf::Contig';

  # check for defined mismatch
  return 0 if((defined $contig->name      ? 1 : 0) != (defined $contig1->name     ? 1 : 0));
  return 0 if((defined $contig->length    ? 1 : 0) != (defined $contig1->length   ? 1 : 0));

  return 0 if defined $contig1->name && $contig->name ne $contig1->name;
  return 0 if defined $contig1->length && $contig->length != $contig1->length;

  return 1 if(!defined $contig->checksum || !defined $contig1->checksum);
  return 0 if $contig->checksum ne $contig1->checksum;
  return 1;
}


1;
