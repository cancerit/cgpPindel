#!/usr/bin/perl
# Copyright (c) 2014-2022 Genome Research Ltd
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


BEGIN {
  use Cwd qw(abs_path);
  use File::Basename;
  unshift (@INC,dirname(abs_path($0)).'/../lib');
};

use strict;
use warnings FATAL => 'all';
use autodie qw(:all);
use Getopt::Long;
use Pod::Usage qw(pod2usage);
use Const::Fast qw(const);
use File::Which qw(which);
use Vcf;
use Sanger::CGP::Pindel;
use PCAP::Cli;


const my $PINDEL_RAW_MIN => 3;
const my $MIN_CALL_FRAC => 0.05;
const my $GFF_HEAD => "##gff-version 3\n##cgpPindel-version: %s\n";
const my $GFF_FORMAT => "%s\tpindel_np_from_vcf\tindel\t%d\t%d\t.\t.\tSAMPLE_COUNT=%d,%s\n";
const my $BED_HEAD => "##cgpPindel-version: %s\n";
const my $BED_FORMAT => "%s\t%d\t%d\t%s\tSAMPLE_COUNT=%s,%s\n";

use Data::Dumper;

my $options = &get_options;
my $calls = collate_calls($options);
my $out_file = output($options, $calls);
bgzip_tabix($out_file, $options->{'r'});

sub collate_calls {
  my $options = shift;
  my $all_calls = {};
  for my $file (@{$options->{'files'}}) {
    my ($new_sample, $new_calls) = parse_vcf($options, $file);
    condense_calls($all_calls, $new_sample, $new_calls, $options->{'r'});
  }
  return $all_calls;
}

sub condense_calls {
  my ($all_calls, $new_sample, $new_calls, $is_range) = @_;
  my $new;
  if($is_range == 1) {
    $new = bed_condense_calls($all_calls, $new_sample, $new_calls);
  }
  else {
    $new = gff_condense_calls($all_calls, $new_sample, $new_calls);
  }
  warn "New rows: $new\t($new_sample)\n";
  return 1;
}

sub bed_condense_calls {
  my ($all_calls, $new_sample, $new_calls) = @_;
  my $new = 0;

  for my $chr(keys %{$new_calls}) {
    for my $rs(keys %{$new_calls->{$chr}}) {
      for my $re(keys %{$new_calls->{$chr}->{$rs}}) {
        for my $pc(keys %{$new_calls->{$chr}->{$rs}->{$re}}) {
          my ($final_calls, $final_depth) = (0,0);
          for my $counts(@{$new_calls->{$chr}->{$rs}->{$re}->{$pc}}) {
            my ($uniq_calls, $uniq_depth) = @{$counts};
            $final_calls += $uniq_calls; # as pindel doesn't share reads between events
            # just take the highest as mixed events possible within poly/rep range
            $final_depth = $uniq_depth if($uniq_depth > $final_depth);
          }
          $new++ unless(exists $all_calls->{$chr}->{$rs}->{$re}->{$pc});
          push @{$all_calls->{$chr}->{$rs}->{$re}->{$pc}}, "$new_sample=$final_calls|$final_depth";
        }
      }
    }
  }
  return $new;
}

sub gff_condense_calls {
  my ($all_calls, $new_sample, $new_calls) = @_;
  my $new = 0;
  for my $chr(keys %{$new_calls}) {
    for my $pos(keys %{$new_calls->{$chr}}) {
      my ($final_calls, $final_depth) = (0,0);
      for my $counts(@{$new_calls->{$chr}->{$pos}}) {
        my ($uniq_calls, $uniq_depth) = @{$counts};
        $final_calls += $uniq_calls; # as pindel doesn't share reads between events (in this version)
        $final_depth = $uniq_depth if($uniq_depth > $final_depth); # as the largest event starting here will have the highest depth
      }
      $new++ unless(exists $all_calls->{$chr}->{$pos});
      push @{$all_calls->{$chr}->{$pos}}, "$new_sample=$final_calls|$final_depth";
    }
  }
  return $new;
}

sub parse_vcf {
  my ($options, $file) = @_;
  my $calls = {};
  my $samp_id = $options->{'samp_id'};

  my $vcf = Vcf->new(file=>$file);
  $vcf->parse_header();
  my $req_sample = ($vcf->get_header_line(key=>'SAMPLE', ID=>$samp_id))[0][0]{'SampleName'};
  $vcf->set_samples(include=>[$samp_id]);
  my ($total, $kept, $new) = (0,0, 0);
  while(my $d = $vcf->next_data_hash) {
    $total++;
    next if($d->{'REF'} eq $d->{'ALT'});
    my $gtypes = $d->{'gtypes'}->{$samp_id};

    next if(($gtypes->{'PP'} + $gtypes->{'NP'}) < $PINDEL_RAW_MIN);
    my $uniq_calls = $gtypes->{'PU'} + $gtypes->{'NU'};
    my $uniq_depth = $gtypes->{'PR'} + $gtypes->{'NR'};
    $uniq_depth ||= 1;
    next if($uniq_calls / $uniq_depth < $MIN_CALL_FRAC);
    if($options->{'r'} == 1) {
      push @{$calls->{$d->{'CHROM'}}->{$d->{'INFO'}->{'RS'}}->{$d->{'INFO'}->{'RE'}}->{$d->{'INFO'}->{'PC'}}}, [$uniq_calls,$uniq_depth];
    }
    else {
      my $gff_corrected_pos = $d->{'POS'} + 1;
      push @{$calls->{$d->{'CHROM'}}->{$gff_corrected_pos}}, [$uniq_calls,$uniq_depth];
    }
    $kept++;
  }
  warn "Retained: $kept/$total\n";
  return ($req_sample, $calls);
}

sub output {
  my ($options, $calls_in) = @_;
  my $ret_file;
  if ($options->{'r'} == 1) {
    $ret_file = output_bed($options, $calls_in);
  }
  else {
    $ret_file = output_gff($options, $calls_in);
  }
  return $ret_file;
}

sub output_bed {
  my ($options, $calls_in) = @_;
  my %calls = %{$calls_in}; # mainly for readability
  my $raw_bed = "$options->{output}.bed";
  open my $bed, '>', $raw_bed;
  print $bed sprintf $BED_HEAD, Sanger::CGP::Pindel->VERSION or die "Failed to write to $raw_bed";
  print $bed "##$options->{cmd}\n" or die "Failed to write to $raw_bed";
  for my $chr(sort keys %calls) {
    for my $rs(sort {$a<=>$b} keys %{$calls{$chr}}) {
      for my $re(sort {$a<=>$b} keys %{$calls{$chr}{$rs}}) {
        for my $pc(sort keys %{$calls{$chr}{$rs}{$re}}) {
          my @sample_sets = sort @{$calls{$chr}{$rs}{$re}{$pc}};
          print $bed sprintf $BED_FORMAT, $chr, $rs-1, $re, $pc, scalar @sample_sets, (join q{,}, @sample_sets) or die "Failed to write to $raw_bed";
        }
      }
    }
  }
  close $bed;
  return $raw_bed;
}

sub output_gff {
  my ($options, $calls_in) = @_;
  my %calls = %{$calls_in}; # mainly for readability

  my $raw_gff3 = "$options->{output}.gff3";
  open my $gff3, '>', $raw_gff3;
  print $gff3 sprintf $GFF_HEAD, Sanger::CGP::Pindel->VERSION or die "Failed to write to $raw_gff3";
  print $gff3 "##$options->{cmd}\n" or die "Failed to write to $raw_gff3";
  for my $chr(sort keys %calls) {
    for my $pos(sort {$a<=>$b} keys %{$calls{$chr}}) {
      my @sample_sets = sort @{$calls{$chr}{$pos}};
      print $gff3 sprintf $GFF_FORMAT, $chr, $pos, $pos, scalar @sample_sets, (join q{,}, @sample_sets) or die "Failed to write to $raw_gff3";
    }
  }
  close $gff3;
  return $raw_gff3;
}

sub bgzip_tabix {
  my ($out_file, $is_range) = @_;

  my $out_gz = "$out_file.gz";
  my $tbi = "$out_gz.tbi";
  unlink $out_gz if(-e $out_gz);
  unlink $tbi if(-e $tbi);

  # autodie makes this nice and clean
  my $bgzip = which('bgzip');
  system($bgzip, $out_file);

  my $tabix = which('tabix');
  my $mode = $is_range == 1 ? 'bed' : 'gff';
  system($tabix, '-p', $mode, $out_gz);
  return 1;
}

sub get_options {
  my %opts = ('r' => 0);
  $opts{'cmd'} = join " ", $0, @ARGV;
  GetOptions( 'h|help'      => \$opts{'h'},
              'm|man'       => \$opts{'m'},
              'v|version'   => \$opts{'v'},
              'r|range'     => \$opts{'r'},
              'o|output=s'  => \$opts{'output'},
              's|samp_id=s' => \$opts{'samp_id'},
  ) or pod2usage(2);

  pod2usage(-verbose => 1) if(defined $opts{'h'});
  pod2usage(-verbose => 2) if(defined $opts{'m'});

  if(defined $opts{'v'}){
    print 'Version: ',Sanger::CGP::Pindel->VERSION,"\n";
    exit;
  }

  pod2usage(-message => q{ERROR: 'output' must be defined.}, -verbose => 1) unless(defined $opts{'output'});
  pod2usage(-message => q{ERROR: 'samp_id' must be defined.}, -verbose => 1) unless(defined $opts{'samp_id'});

  $opts{'files'} = \@ARGV;

  for my $file(@{$opts{'files'}}) {
    PCAP::Cli::file_for_reading($file, $file);
  }

  return \%opts;
}

##FORMAT=<ID=PR,Number=1,Type=Integer,Description="Total mapped reads on the positive strand">
##FORMAT=<ID=NR,Number=1,Type=Integer,Description="Total mapped reads on the negative strand">
##FORMAT=<ID=PU,Number=1,Type=Integer,Description="Unique calls on the positive strand">
##FORMAT=<ID=NU,Number=1,Type=Integer,Description="Unique calls on the negative strand">
##SAMPLE_COUNT=2,PD4975b=55|148,PD6722b=69|178

__END__

=head1 pindel_np_from_vcf.pl

Generate a pinel normal panel gff3 from VCF inputs

=head1 SYNOPSIS

pindel_np_from_vcf.pl [options] *.vcf[.gz]

  Required options:

    -output   -o  File stub to write final result to
                   - appropriate extention will be appended.
    -samp_id  -s  Which sample ID should be used
                   - output from pindel.pl tags all data
                     as TUMOUR or NORMAL
    -range    -r  Generate range based normal panel
                   - without a gff3.gz file is generated (change start only, original)
                   - with a bed.gz file is generated (range start/end based)

  Other:
    -help     -h  Brief help message.
    -man      -m  Full documentation.
    -version  -v  Version
