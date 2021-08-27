#!/usr/bin/perl
# Copyright (c) 2014-2021 Genome Research Ltd
#
# Author: CASM/Cancer IT <cgphelp@sanger.ac.uk>
#
# This file is part of cgpPindel.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
# 1. The usage of a range of years within a copyright statement contained within
# this distribution should be interpreted as being equivalent to a list of years
# including the first and last year specified and all consecutive years between
# them. For example, a copyright statement that reads ‘Copyright (c) 2005, 2007-
# 2009, 2011-2012’ should be interpreted as being identical to a statement that
# reads ‘Copyright (c) 2005, 2007, 2008, 2009, 2011, 2012’ and a copyright
# statement that reads ‘Copyright (c) 2005-2012’ should be interpreted as being
# identical to a statement that reads ‘Copyright (c) 2005, 2006, 2007, 2008,
# 2009, 2010, 2011, 2012’.
#
use strict;
use warnings FATAL => 'all';
use autodie qw(:all);
use Cwd qw(abs_path);
use Pod::Usage qw(pod2usage);
use FindBin qw($Bin);
use lib "$Bin/../lib";

use File::Basename;
use Getopt::Long;
use Data::UUID;

use Sanger::CGP::Pindel::Implement;
use Sanger::CGP::Pindel::OutputGen::VcfCohortConverter;
use Sanger::CGP::Vcf::Contig;
use Sanger::CGP::Vcf::Sample;
use Sanger::CGP::Pindel::OutputGen::PindelRecordParser;
use PCAP::Cli;
use PCAP::Bam::Bas;

use Bio::DB::HTS;

use Data::Dumper;

{
  my $options = setup();
  my $contigs = contig_setup($options);
  my ($hts_by_sample, $bas_by_sample, $vcfsamp_by_sample) = hts_and_sample_parse($options);

  my @samples = sort keys %{$hts_by_sample};
  my $record_converter = Sanger::CGP::Pindel::OutputGen::VcfCohortConverter->new(
    -contigs => [values %$contigs],
    -samples => \@samples,
    -hts_set => $hts_by_sample,
    -bas_set => $bas_by_sample,
    -all => $options->{'all'},
    -badloci => $options->{'badloci'},
  );
  my $input_source = basename($0). '_v'. Sanger::CGP::Pindel->VERSION;
  my $header = $record_converter->gen_header($options->{'reference'}, $input_source, $vcfsamp_by_sample, $options);
  my $out_fh = $options->{'output'};
  print $out_fh $header;
  for my $in_file(@{$options->{'input'}}) {
    process_pindel_file($options, $in_file, $out_fh, $record_converter);
  }
  close $options->{'output'} if($options->{'output'} ne \*STDOUT);
}

sub process_pindel_file {
  my ($options, $pindel_file, $out_fh, $record_converter) = @_;
  my $prp = Sanger::CGP::Pindel::OutputGen::PindelRecordParser->new(
    -path => $pindel_file,
    -fai => Bio::DB::HTS::Fai->load($options->{'reference'}),
    -noreads => 1,

  );
  my $uuid_gen = Data::UUID->new;
my $records = 0;
  while(my $record = $prp->next_record) {
    $record->id($uuid_gen->to_string($uuid_gen->create));
#last if($records > 40);
    print $out_fh $record_converter->gen_record($record);
#$records++;
  }
}

sub hts_and_sample_parse {
  my $options = shift;
  my %hts_by_sample;
  my %bas_by_sample;
  my %vcf_samples;
  for my $hts_f(@{$options->{'hts_files'}}) {
    my $bas_f = $hts_f.'.bas';
    die "ERROR: Failed to find colocated bas file at: $bas_f" unless(-e $bas_f);
    my $hts = Bio::DB::HTS->new(-bam => $hts_f, -fasta => $options->{reference});
    my ($sample, $platform);
    foreach my $line (split(/\n/,$hts->header->text)) {
      next unless($line =~ m/^\@RG/);
      chomp $line;
      ($sample) = $line =~ m/SM:([^\t]+)/;
      ($platform) = $line =~ /PL:([^\t]+)/;
      last if(defined $sample);
    }
    $hts_by_sample{$sample} = $hts_f;
    $bas_by_sample{$sample} = PCAP::Bam::Bas->new($bas_f);
    die sprintf "Failed to find sample name in \@RG header lines of %s", $hts->hts_path unless(defined $sample);
    $vcf_samples{$sample} = Sanger::CGP::Vcf::Sample->new(
        -name => $sample,
        -study => $options->{'project'},
        -platform => $platform,
        -seq_protocol => $options->{'protocol'},
        -description => $sample
      );

  }
  return (\%hts_by_sample, \%bas_by_sample, \%vcf_samples);
}

sub _contig_parse {
  my ($hts, $hts_contigs, $fixed_contigs) = @_;
  my ($assembly_out, $species_out);
  foreach my $line (split(/\n/,$hts->header->text)){
    next unless($line =~ /^\@SQ/);
    my ($name) = $line =~ /SN:([^\t]+)/;
    my ($length) = $line =~ /LN:([^\t]+)/;
    my ($assembly) = $line =~ /AS:([^\t]+)/;
    my ($species) = $line =~ /SP:([^\t]+)/;
    if(defined $assembly && !defined $assembly_out) {
      $assembly_out = $assembly;
    }
    if(defined $species && !defined $species_out) {
      $species_out = $species;
    }
    if($fixed_contigs) {
      if(!exists $hts_contigs->{$name}) {
        die sprintf "ERROR: Found contig %s in file %s but this is not consistent with other files.", $name, $hts->hts_path;
      }
    }
    else {
      $hts_contigs->{$name} = $length;
    }
  }
  return ($assembly_out, $species_out);
}

sub contig_setup {
  my $options = shift;
  my @hts; # to store objects
  my $hts_contigs = {}; # to store contig list
  my ($assembly, $species);

  $assembly = $options->{'assembly'} if(defined $options->{'assembly'});
  $species = join q{ }, @{$options->{'species'}} if(defined $options->{'species'} && @{$options->{'species'}} > 0);

  my @paths;
  my $not_first = 0;
  for my $hts_f(@{$options->{'hts_files'}}) {
    my $hts = Bio::DB::HTS->new(-bam => $hts_f, -fasta => $options->{reference});
    my ($assembly_found, $species_found) = _contig_parse($hts, $hts_contigs, $not_first++);
    if(! defined $assembly && defined $assembly_found) {
      $assembly = $assembly_found;
    }
    if(defined $assembly && defined $assembly_found && $assembly ne $assembly_found) {
      if(defined $options->{'assembly'}){
        die sprintf "Assembly defined on command line (%s) doesn't match that found in %s (%s)", $options->{'assembly'}, $hts->hts_path, $assembly_found;
      }
      else {
        die sprintf "Assembly defined in %s (%s) doesn't match that found in previous files (%s):%s\n", $hts->hts_path, $assembly_found, $assembly, join(qq{\n\t}, @paths);
      }
    }
    if(! defined $species && defined $species_found) {
      $species = $species_found;
    }
    if(defined $species && defined $species_found && $species ne $species_found) {
      if(defined $options->{'species'}){
        die sprintf "Species defined on command line (%s) doesn't match that found in %s (%s)", $options->{'species'}, $hts->hts_path, $species_found;
      }
      else {
        die sprintf "Species defined in %s (%s) doesn't match that found in previous files (%s):%s\n", $hts->hts_path, $species_found, $species, join(qq{\n\t}, @paths);
      }
    }
    push @paths, $hts->hts_path;
  }
  die "No assembly defined in BAM/CRAM headers please specify in command options." unless(defined $assembly);
  die "No species defined in BAM/CRAM headers please specify in command options." unless(defined $species);
  for my $name(keys %{$hts_contigs}) {
    my $contig = Sanger::CGP::Vcf::Contig->new(
      -name => $name,
      -length => $hts_contigs->{$name},
      -assembly => $assembly,
      -species => $species
    );
    $hts_contigs->{$name} = $contig;
  }
  return $hts_contigs;
}

sub setup{
  my %opts = (
    'cmd' => join(" ", $0, @ARGV),
    'all' => 0,
  );
  GetOptions( 'h|help' => \$opts{'h'},
              'm|man' => \$opts{'m'},
              'v|version' => \$opts{'v'},
              'o|output:s' => \$opts{'output'},
              'r|ref=s' => \$opts{'reference'},
              'i|input=s@' => \$opts{'input'},
              'prj|project:s' => \$opts{'project'},
              'p|prot:s' => \$opts{'protocol'},
              'as|assembly:s' => \$opts{'assembly'},
              'sp|species=s{0,}' => \@{$opts{'species'}},
              'pp|parent:s' => \$opts{'pp'},
              'a|all' => \$opts{'all'},
              'b|badloci:s' => \$opts{'badloci'},
  );

  if(defined $opts{'v'}){
    printf "Version: %s\n", Sanger::CGP::Pindel::Implement->VERSION;
    exit;
  }

  pod2usage(-verbose => 1) if(defined $opts{'h'});
  pod2usage(-verbose => 2) if(defined $opts{'m'});

  PCAP::Cli::file_for_reading('ref', $opts{'reference'});
  my @full_inputs;
  for my $if(@{$opts{'input'}}) {
    if($if =~ s/%/*/g) {
      for my $gf(glob $if) {
        push @full_inputs, $gf if(-s $gf);
      }
    }
    else {
      push @full_inputs, $if;
    }
  }
  for my $if(@full_inputs) {
    PCAP::Cli::file_for_reading('input (possibly expanded)', $if);
  }
  $opts{'input'} = \@full_inputs;

  if($opts{'output'}) {
    open my $fh, '>', $opts{'output'};
    $opts{'output'} = $fh;
  }
  else { $opts{'output'} = \*STDOUT; }

  # add hts_files from the remains of @ARGV
  Sanger::CGP::Pindel::Implement::cohort_files(\%opts);

  return \%opts;
}


__END__

=head1 NAME

pindelCohort_to_vcf.pl - Takes raw Pindel files and a set of bam files to produces a collated vcf file.

=head1 SYNOPSIS

pindelCohort_to_vcf.pl [options] SAMPLE1.bam [SAMPLE2.bam ...]

  Required parameters:
    -ref       -r   File path to the reference file used to provide the coordinate system.
    -input     -i   Files to read in, repeatable or '%' wildcard

  Optional parameters:
    -output    -o   File path to output to. Defaults to STDOUT.
    -all       -a   Generate VAF for all samples, even when not seen by Pindel.
    -badloci   -b   Tabix indexed BED file of locations reject as events
                     - e.g. hi-seq depth from UCSC
    -project   -prj String representing the project data is from.
    -prot      -p   String representing the sequencing protocol (e.g. genomic, targeted, RNA-seq).
    -assembly  -as  Reference assembly name, used when not found in BAM/CRAM headers.
    -species   -sp  Species name, used when not found in BAM/CRAM headers.
    -parent    -pp  Process information from parent program (where one exists)

  Other:
    -help      -h   Brief help message.
    -man       -m   Full documentation.
    -version   -v   Prints the version number.

=head1 OPTIONS

=over 8

=item B<-input>

File path(s) to read. Accepts only raw pindel files, repeat or '%' as wildcard for multiple.

=item B<-output>

File path to output data. If this option is omitted the script will attempt to write to STDOUT.


=item B<-project>

String identifying the project to which the samples belong.

=item B<-prot>

String representing the sequencing protocol (e.g. genomic, targeted, RNA-seq).

=item B<-assembly>

Reference assembly name, used when not found in BAM headers.  Validated against header if both are present.

=item B<-species>

Species name, used when not found in BAM headers.  Validated against header if both are present.

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=item B<-version>

Prints the version number and exits.

=back

=head1 DESCRIPTION

B<pindelCohort_to_vcf.pl> will attempt to generate a vcf file from a set of Pindel output files.

=cut
