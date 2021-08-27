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
BEGIN {
  use Cwd qw(abs_path cwd);
  use File::Basename;
  unshift (@INC,dirname(abs_path($0)).'/../lib');
};

use strict;
use warnings FATAL => 'all';
use autodie qw(:all);
use Const::Fast qw(const);
use File::Copy;
use File::Path qw(remove_tree make_path);

use PCAP::Cli;
use Sanger::CGP::Pindel::Implement;

const my %INDEX_MAX => (
                  'input'   => 1,
                  'pindel'  => -1,
                  'parse' => 1, # reads all pout, makes raw-vcf and splits to even sized for blat
                  'blat' => -1,
                  'concat' => 1,
                  );
const my @VALID_PROCESS => keys %INDEX_MAX;

{
  my $options = setup();
  my $threads = PCAP::Threaded->new($options->{'threads'});

  # start processes here (in correct order obviously), add conditions for skipping based on 'process' option
  if(!exists $options->{'process'} || $options->{'process'} eq 'input') {
    Sanger::CGP::Pindel::Implement::input_cohort($options)
  }
  if(!exists $options->{'process'} || $options->{'process'} eq 'pindel') {
    my $jobs = Sanger::CGP::Pindel::Implement::determine_jobs($options); # method still needed to populate info
    $jobs = $options->{'limit'} if(exists $options->{'limit'} && defined $options->{'limit'});
    $threads->add_function('pindel', \&Sanger::CGP::Pindel::Implement::pindel);
    $threads->run($jobs, 'pindel', $options);
  }
  if(!exists $options->{'process'} || $options->{'process'} eq 'parse') {
    Sanger::CGP::Pindel::Implement::parse($options);
  }
  $options->{'split_files'} = Sanger::CGP::Pindel::Implement::split_files($options) unless(exists $options->{'split_files'});
  if(!exists $options->{'process'} || $options->{'process'} eq 'blat') {
    my $jobs = scalar @{$options->{'split_files'}};
    $jobs = $options->{'limit'} if(exists $options->{'limit'} && defined $options->{'limit'});
    $threads->add_function('blat', \&Sanger::CGP::Pindel::Implement::blat);
    $threads->run($jobs, 'blat', $options);
  }
  if(!exists $options->{'process'} || $options->{'process'} eq 'concat') {
    Sanger::CGP::Pindel::Implement::concat($options);
    cleanup($options) unless($options->{'debug'});
  }
}

sub cleanup {
  my $options = shift;
  my $tmpdir = $options->{'tmp'};
  move(File::Spec->catdir($tmpdir, 'logs'), File::Spec->catdir($options->{'outdir'}, 'logs')) || die $!;
  remove_tree $tmpdir if(-e $tmpdir);
	return 0;
}

sub index_check {
  my $opts = shift;
  my $max_files = @{$opts->{'hts_files'}};
  if(exists $opts->{'process'}) {
    PCAP::Cli::valid_process('process', $opts->{'process'}, \@VALID_PROCESS);
    if(exists $opts->{'index'}) {
      my @valid_seqs = Sanger::CGP::Pindel::Implement::valid_seqs($opts);
      my $refs = scalar @valid_seqs;

      my $max = $INDEX_MAX{$opts->{'process'}};
      if($max == -1){
        if(exists $opts->{'limit'}) {
          $max = $opts->{'limit'} > $refs ? $refs : $opts->{'limit'};
        } else {
          if($opts->{'process'} eq 'input') {
            $max = $max_files;
          }
          elsif($opts->{'process'} eq 'blat') {
            $opts->{'split_files'} = Sanger::CGP::Pindel::Implement::split_files($opts);
            $max = scalar @{$opts->{'split_files'}};
          }
          else {
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
  $opts->{pad} = 1 unless(exists $opts->{pad} && defined $opts->{pad});

  # add hts_files from the remains of @ARGV
  Sanger::CGP::Pindel::Implement::cohort_files($opts);
  index_check($opts);

  return $opts;
}

__END__

=head1 pindelCohort.pl

Similar to pindel.pl but processes 1 sample. References to BAM can be replaced with CRAM.

=head1 SYNOPSIS

pindelCohort.pl [options] sample1.bam

  Required parameters:
    -outdir    -o   Folder to output result to.
    -reference -r   Path to reference genome file *.fa[.gz]

  Optional
    -pad            Multiples (>=1) of max readlength to pad blat target seq with [default 1]
    -seqtype   -st  Sequencing protocol, expect all input to match [WGS]
    -assembly  -as  Name of assembly in use
                     -  when not available in BAM header SQ line.
    -species   -sp  Species
                     -  when not available in BAM header SQ line.
    -exclude   -e   Exclude this list of ref sequences from processing, wildcard '%'
                     - comma separated, e.g. NC_007605,hs37d5,GL%
    -badloci   -b   Tabix indexed BED file of locations to not accept as anchors or valid events
                     - e.g. hi-seq depth from UCSC
    -cpus      -c   Number of cores to use. [1]
                     - recommend max 4 during 'input' process.
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

  input
  pindel - index available
  parse
  blat - index available
  concat

=item B<-index>

Possible index ranges for processes above are:

  ?

If you want STDOUT/ERR to screen ensure index is set even for single job steps.

=back
