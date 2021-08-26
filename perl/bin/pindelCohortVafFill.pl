#!/usr/bin/env perl
# Copyright (c) 2014-2021
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
use Getopt::Long;
use File::Spec::Functions;
use File::Path qw(make_path remove_tree);
use Const::Fast qw(const);
use File::Copy qw(move);

use PCAP::Cli;
use PCAP::Threaded;
use Sanger::CGP::Pindel::Implement;

use Data::Dumper;

const my %INDEX_MAX => (
                  'split'   => 1,
                  'fill'  => -1, # depends on number of splits
                  'bams' => -1, # may want to make multi but probably not necessary
                  'finalise' => 1, # merging of split files
                  );
const my @VALID_PROCESS => keys %INDEX_MAX;

{
  my $options = setup();
  my $threads = PCAP::Threaded->new($options->{threads});

  if(!exists $options->{process} || $options->{process} eq 'split') {
    Sanger::CGP::Pindel::Implement::cohort_split($options);
    vaf_fill_seqdata($options);
  }

  if(!exists $options->{process} || $options->{process} eq 'fill') {
    for my $f(glob(catfile($options->{'split_dir'}, '*.vcf.gz'))) {
      push @{$options->{split_files}}, $f if($f =~ m{/\d+.vcf.gz$});
    }
    $threads->add_function('fill', \&Sanger::CGP::Pindel::Implement::fill_split_vaf);
    $threads->run(scalar @{$options->{split_files}}, 'fill', $options);
  }

  if(!exists $options->{process} || $options->{process} eq 'bams') {
    $threads->add_function('bams', \&Sanger::CGP::Pindel::Implement::merge_vaf_bams);
    $threads->run(scalar @{$options->{primary_hts}}, 'bams', $options);
  }

  if(!exists $options->{process} || $options->{process} eq 'finalise') {
    Sanger::CGP::Pindel::Implement::fill_vcf_merge($options);
    if(!$options->{debug}) {
      move(catdir($options->{tmp}, 'logs'), catdir($options->{output}, 'logs'));
      remove_tree($options->{tmp});
    }
  }
}

sub vaf_fill_seqdata {
  my ($options) = @_;
  my $bwa_files = $options->{bwa_file_list};
  return if (-e $bwa_files);
  open my $FH, '>', $bwa_files;
  for my $f(@{$options->{primary_hts}}) {
    print $FH qq{$f\n};
  }
  close $FH;
}

sub setup {
  my %opts = (
    'size' => 10000,
    'name' => 'cohort',
    'primary_hts' => [],
    'secondary_hts' => [],
    'threads' => 1,
  );
  GetOptions( 'h|help' => \$opts{h},
              'm|man' => \$opts{'m'},
              'v|version' => \$opts{v},
              'f|input=s' => \$opts{input},
              'r|ref=s' => \$opts{ref},
              'o|output=s' => \$opts{output},
              's|size:i' => \$opts{size},
              'n|name:s' => \$opts{name},
              'c|cpus:i' => \$opts{threads},
              'p|process:s' => \$opts{process},
              'i|index:i' => \$opts{index},
              'l|limit:i' => \$opts{limit},
              'a|abort' => \$opts{abort},
              'd|data=s' => \$opts{data},
              'debug' => \$opts{debug}
  );

  if(defined $opts{v}) {
    printf "Version: %s\n", Sanger::CGP::Pindel::Implement->VERSION;
    exit;
  }
  pod2usage(-verbose => 1) if(defined $opts{h});
  pod2usage(-verbose => 2) if(defined $opts{'m'});


  delete $opts{'process'} unless(defined $opts{'process'});
  delete $opts{'index'} unless(defined $opts{'index'});
  delete $opts{'limit'} unless(defined $opts{'limit'});

  for my $param(qw(input ref output data)) {
    pod2usage(-verbose => 1, -message=> sprintf('ERROR: -%s must be defined', $param), -exit => 2) unless(defined $opts{$param});
  }

  my $final_logs = catdir($opts{output}, 'logs');
  if(-e $final_logs) {
    warn "WARN: $final_logs directory exists suggesting completed analysis.\n";
    if($opts{abort}) {
      die "EXIT as -abort in effect.\n";
    }
    else {
      warn "WARN: Continuing, set -abort to prevent reprocessing of completed data.\n";
    }
  }

  PCAP::Cli::file_for_reading('input', $opts{input});
  PCAP::Cli::file_for_reading('data', $opts{data});
  PCAP::Cli::file_for_reading('ref', $opts{ref});

  if(@ARGV) {
    die "ERROR: No positional arguments expected."
  }

  open my $D, '<', $opts{data};
  while(my $l = <$D>) {
    chomp $l;
    my ($primary, $secondary);
    if($l =~ m/^([^\t]+)\t([^\t]+)$/) {
      ($primary, $secondary) = ($1, $2);
    }
    else {
      die "ERROR: file '$opts{data}' is malformed\n";
    }
    PCAP::Cli::file_for_reading('bam/cram files', $primary);
    PCAP::Cli::file_for_reading('bam/cram files', $secondary);
    push @{$opts{primary_hts}}, $primary;
    push @{$opts{secondary_hts}}, $secondary;
  }
  close $D;
  my $sample_count = @{$opts{primary_hts}};
  if($opts{size} < $sample_count * 5) {
    $opts{size} = $sample_count * 5;
    warn sprintf "WARNING: -size has been automatically increased to %d (5x sample number), see '-help'\n", $opts{size};
  }

  $opts{tmp} = catdir($opts{output}, 'tmpCohortVafFill');
  make_path($opts{tmp});

  $opts{split_dir} = catdir($opts{tmp}, 'split');
  make_path($opts{split_dir});
  $opts{fill_dir} = catdir($opts{tmp}, 'fill');
  make_path($opts{fill_dir});
  $opts{bwa_file_list} = catfile($opts{tmp}, 'bwa_files.lst');

  make_path(catdir($opts{tmp}, 'logs'));

  return \%opts;
}

__END__

=head1 NAME

pindelCohortVafFill.pl - Takes merged cohort VCF and fills in gaps in farm friendly manner.

=head1 SYNOPSIS

pindelCohortVafFill.pl [options] -i ... -o ... -r ... -d ...

  Required parameters:
    -input     -f   VCF file to read in.
    -output    -o   Workspace directory and final output.
    -ref       -r   File path to the reference file used to provide the coordinate system.
    -data      -d   File containing list of sequence data files for all samples used in "-input"
                    - format: tab separated BWA mapping followed by pindel_cohort reads, one sample per line.

                        sample_A_bwa.bam<TAB>sample_A_pindel.bam
                        sample_B_bwa.bam<TAB>sample_B_pindel.bam

  Optional parameters:
    -name      -n   Stub name for final output files [$output/cohort...]
    -size      -s   Number of Sample/event combinations per-file when processing [10000]
                     - Automatically increased to a minimum of 5x number of samples (for efficiency).
    -abort     -a   Abort noisily if data appears to have been processed (silent exit otherwise)
    -cpus      -c   Number of cores to use. [1]


  Targeted processing (further detail under OPTIONS):
    -process   -p   Only process this step then exit, optionally set -index
    -index     -i   Optionally restrict '-p' to single job
    -limit     -l   Use '-p fill -i N -l M' - do not declare '-c'
                      M: maximum number of scattered processes on multiple hosts
                      N: This instance, 1..M

  Other:
    -help      -h   Brief help message.
    -man       -m   Full documentation.
    -version   -v   Prints the version number.
    -debug          Don't cleanup anything

=head1 DESCRIPTION

B<pindelCohortVafFill.pl> Takes merged cohort VCF and fills in gaps in farm friendly manner.

One additional file is generated containg events where no fill-in is required.

Each original sample BAM/CRAM (normally BWA mapping) needs to be paired with the one generated by pindelCohort.pl

Each sample BAM/CRAM also needs colocated index and bas file.

=cut
