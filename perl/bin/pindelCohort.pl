#!/usr/bin/perl

########## LICENCE ##########
# Copyright (c) 2020 Genome Research Ltd.
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

BEGIN {
  use Cwd qw(abs_path cwd);
  use File::Basename;
  unshift (@INC,dirname(abs_path($0)).'/../lib');
};

use strict;
use warnings FATAL => 'all';
use autodie qw(:all);
use Const::Fast qw(const);

use PCAP::Cli;
use Sanger::CGP::Pindel::Implement;

use Data::Dumper;

const my @VALID_PROCESS => qw(input pindel ???);
my %index_max = ( 'input'   => -1,
                  'pindel'  => -1,
                  '???' => 1,
                  );

{
  my $options = setup();
  my $threads = PCAP::Threaded->new($options->{'threads'});
  &PCAP::Threaded::disable_out_err if(exists $options->{'index'});

  # register any process that can run in parallel here
  $threads->add_function('input', \&Sanger::CGP::Pindel::Implement::input_cohort);
  $threads->add_function('pindel', \&Sanger::CGP::Pindel::Implement::pindel);

  # start processes here (in correct order obviously), add conditions for skipping based on 'process' option
  if(!exists $options->{'process'} || $options->{'process'} eq 'input') {
    my $jobs = scalar @{$options->{'hts_files'}};
    $jobs = $options->{'limit'} if(exists $options->{'limit'} && defined $options->{'limit'});
    $threads->run($jobs, 'input', $options);
  }
  if(!exists $options->{'process'} || $options->{'process'} eq 'pindel') {
    my $jobs = Sanger::CGP::Pindel::Implement::determine_jobs($options); # method still needed to populate info
    $jobs = $options->{'limit'} if(exists $options->{'limit'} && defined $options->{'limit'});
    $threads->run($jobs, 'pindel', $options);
  }
}

sub index_check {
  my $opts = shift;
  my $max_files = @{$opts->{'hts_files'}};
  if(exists $opts->{'process'}) {
    PCAP::Cli::valid_process('process', $opts->{'process'}, \@VALID_PROCESS);
    if(exists $opts->{'index'}) {
      my @valid_seqs = Sanger::CGP::Pindel::Implement::valid_seqs($opts);
      my $refs = scalar @valid_seqs;

      my $max = $index_max{$opts->{'process'}};
      if($max==-1){
        if(exists $opts->{'limit'}) {
          $max = $opts->{'limit'} > $refs ? $refs : $opts->{'limit'};
        } else {
          if($opts->{'process'} eq 'input') {
            $max = $max_files;
          } else {
            $max = $refs;
          }
        }
      }
      if($opts->{'index'} < 1 || $opts->{'index'} > $max) {
        if($opts->{'process'} eq 'input') {
          die "ERROR: based on number of inputs option -index must be between 1 and $max_files\n";
        } else {
          die "ERROR: based on reference and exclude option -index must be between 1 and $refs\n";
        }
      }
      PCAP::Cli::opt_requires_opts('index', $opts, ['process']);
      die "No max has been defined for this process type\n" if($max == 0);
    }
  }
  elsif(exists $opts->{'index'}) {
    die "ERROR: -index cannot be defined without -process\n";
  }
}

sub setup {
  my $opts = Sanger::CGP::Pindel::Implement::shared_setup([],{});

  # add hts_files from the remains of @ARGV
  Sanger::CGP::Pindel::Implement::cohort_files($opts);
  index_check($opts);

  return $opts;
}

__END__

=head1 pindelCohort.pl

Similar to pindel.pl but processes 1 or more samples. References to BAM can ber replaced with CRAM.

=head1 SYNOPSIS

pindelCohort.pl [options] sample1.bam [sample2.bam sample3.bam...]

  Required parameters:
    -outdir    -o   Folder to output result to.
    -reference -r   Path to reference genome file *.fa[.gz]
    -simrep    -s   Full path to tabix indexed simple/satellite repeats.
    -filter    -f   VCF filter rules file (see FlagVcf.pl for details)
    -genes     -g   Full path to tabix indexed coding gene footprints.
    -unmatched -u   Full path to tabix indexed gff3 of unmatched normal panel
                      - see pindel_np_from_vcf.pl

  Optional
    -seqtype   -st  Sequencing protocol, expect all input to match [WGS]
    -assembly  -as  Name of assembly in use
                     -  when not available in BAM header SQ line.
    -species   -sp  Species
                     -  when not available in BAM header SQ line.
    -exclude   -e   Exclude this list of ref sequences from processing, wildcard '%'
                     - comma separated, e.g. NC_007605,hs37d5,GL%
    -badloci   -b   Tabix indexed BED file of locations to not accept as anchors
                     - e.g. hi-seq depth from UCSC
    -skipgerm  -sg  Don't output events with more evidence in normal BAM.
    -cpus      -c   Number of cores to use. [1]
                     - recommend max 4 during 'input' process.
    -softfil   -sf  VCF filter rules to be indicated in INFO field as soft flags
    -limit     -l   When defined with '-cpus' internally thread concurrent processes.
                     - requires '-p', specifically for pindel/pin2vcf steps
    -debug     -d   Don't cleanup workarea on completion.
    -apid      -a   Analysis process ID (numeric) - for cgpAnalysisProc header info
                     - not necessary for external use

  Targeted processing (further detail under OPTIONS):
    -process   -p   Only process this step then exit, optionally set -index
    -index     -i   Optionally restrict '-p' to single job

  Other:
    -help      -h   Brief help message.
    -man       -m   Full documentation.
    -version   -v   Version

  File list can be full file names or wildcard, e.g.

    pindelCohort.pl -c 4 -r some/genome.fa[.gz] -o myout sample1.bam sample2.bam
    or
    pindelCohort.pl -c 4 -r some/genome.fa[.gz] -o myout sample*.bam
    or
    pindelCohort.pl -c 4 -r some/genome.fa[.gz] -o myout sample*.cram

  Please note that colocated index and *.bas files are required.

=head1 OPTIONS

=over 2

=item B<-process>

Available processes for this tool are:

  ???

=item B<-index>

Possible index ranges for processes above are:

  ?

If you want STDOUT/ERR to screen ensure index is set even for single job steps.

=back

