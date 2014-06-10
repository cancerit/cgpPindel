#!/usr/bin/perl

########## LICENCE ##########
# Copyright (c) 2014 Genome Research Ltd.
#
# Author: Keiran Raine <cgpit@sanger.ac.uk>
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
  use Cwd qw(abs_path);
  use File::Basename;
  unshift (@INC,dirname(abs_path($0)).'/../lib');
};

use strict;
use warnings FATAL => 'all';
use autodie qw(:all);

use File::Path qw(remove_tree make_path);
use Getopt::Long;
use File::Spec;
use Pod::Usage qw(pod2usage);
use List::Util qw(first);
use Const::Fast qw(const);
use File::Copy;

use PCAP::Cli;
use Sanger::CGP::Pindel::Implement;

const my @VALID_PROCESS => qw(input pindel pin2vcf merge flag);
my %index_max = ( 'input'   => 2,
                  'pindel'  => -1,
                  'pin2vcf' => -1,
                  'merge'   => 1,
                  'flag'    => 1,);

{
  my $options = setup();
  Sanger::CGP::Pindel::Implement::prepare($options);
  my $threads = PCAP::Threaded->new($options->{'threads'});
  &PCAP::Threaded::disable_out_err if(exists $options->{'index'});

  # register any process that can run in parallel here
  $threads->add_function('input', \&Sanger::CGP::Pindel::Implement::input, exists $options->{'index'} ? 1 : 2);
  $threads->add_function('pindel', \&Sanger::CGP::Pindel::Implement::pindel);
  $threads->add_function('pin2vcf', \&Sanger::CGP::Pindel::Implement::pindel_to_vcf);

  # start processes here (in correct order obviously), add conditions for skipping based on 'process' option
  $threads->run(2, 'input', $options) if(!exists $options->{'process'} || $options->{'process'} eq 'input');

  # count the valid input files, gives constant job count for downstream
  my $jobs = Sanger::CGP::Pindel::Implement::determine_jobs($options) if(!exists $options->{'process'} || first { $options->{'process'} eq $_ } ('pindel', 'pin2vcf'));

  $threads->run($jobs, 'pindel', $options) if(!exists $options->{'process'} || $options->{'process'} eq 'pindel');
  $threads->run($jobs, 'pin2vcf', $options) if(!exists $options->{'process'} || $options->{'process'} eq 'pin2vcf');

  Sanger::CGP::Pindel::Implement::merge_and_bam($options) if(!exists $options->{'process'} || $options->{'process'} eq 'merge');

  if(!exists $options->{'process'} || $options->{'process'} eq 'flag') {
    Sanger::CGP::Pindel::Implement::flag($options);
    cleanup($options);
  }
}

sub cleanup {
  my $options = shift;
  my $tmpdir = $options->{'tmp'};
  move(File::Spec->catdir($tmpdir, 'logs'), File::Spec->catdir($options->{'outdir'}, 'logs')) || die $!;
  remove_tree $tmpdir if(-e $tmpdir);
	return 0;
}

sub setup {
  my %opts;
  pod2usage(-msg  => "\nERROR: Option must be defined.\n", -verbose => 1,  -output => \*STDERR) if(scalar @ARGV == 0);
  $opts{'cmd'} = join " ", $0, @ARGV;
  GetOptions( 'h|help' => \$opts{'h'},
              'm|man' => \$opts{'m'},
              'c|cpus=i' => \$opts{'threads'},
              'r|reference=s' => \$opts{'reference'},
              'o|outdir=s' => \$opts{'outdir'},
              't|tumour=s' => \$opts{'tumour'},
              'n|normal=s' => \$opts{'normal'},
              'e|exclude=s' => \$opts{'exclude'},
              'p|process=s' => \$opts{'process'},
              'i|index=i' => \$opts{'index'},
              'v|version' => \$opts{'version'},
              # these are specifically for pin2vcf
              'sp|species=s' => \$opts{'species'},
              'as|assembly=s' => \$opts{'assembly'},
              'st|seqtype=s' => \$opts{'seqtype'},
              'sg|skipgerm' => \$opts{'skipgerm'},
              # specifically for FlagVCF
              's|simrep=s' => \$opts{'simrep'},
              'f|filters=s' => \$opts{'filters'},
              'g|genes=s' => \$opts{'genes'},
              'u|unmatched=s' => \$opts{'unmatched'},
              'sf|softfil=s' => \$opts{'softfil'},
  ) or pod2usage(2);

  pod2usage(-verbose => 1) if(defined $opts{'h'});
  pod2usage(-verbose => 2) if(defined $opts{'m'});

  if($opts{'version'}) {
    print 'Version: ',Sanger::CGP::Pindel::Implement->VERSION,"\n";
    exit 0;
  }

  PCAP::Cli::file_for_reading('reference', $opts{'reference'});
  PCAP::Cli::file_for_reading('tumour', $opts{'tumour'});
  PCAP::Cli::file_for_reading('normal', $opts{'normal'});
  PCAP::Cli::file_for_reading('simrep', $opts{'simrep'});
  PCAP::Cli::file_for_reading('filters', $opts{'filters'});
  PCAP::Cli::file_for_reading('genes', $opts{'genes'});
  PCAP::Cli::file_for_reading('unmatched', $opts{'unmatched'});
  PCAP::Cli::file_for_reading('softfil', $opts{'softfil'}) if(defined $opts{'softfil'});
  PCAP::Cli::out_dir_check('outdir', $opts{'outdir'});

  delete $opts{'process'} unless(defined $opts{'process'});
  delete $opts{'index'} unless(defined $opts{'index'});

  unless(defined $opts{'exclude'}) {
    delete $opts{'exclude'};
  }

  if(exists $opts{'process'}) {
    PCAP::Cli::valid_process('process', $opts{'process'}, \@VALID_PROCESS);
    if(exists $opts{'index'}) {
      my @valid_seqs = Sanger::CGP::Pindel::Implement::valid_seqs(\%opts);
      my $refs = scalar @valid_seqs;

      my $max = $index_max{$opts{'process'}};
      $max = $refs if($max == -1);

      die "ERROR: based on reference and exclude option index must be between 1 and $refs\n" if($opts{'index'} < 1 || $opts{'index'} > $max);
      PCAP::Cli::opt_requires_opts('index', \%opts, ['process']);

      die "No max has been defined for this process type\n" if($max == 0);

      PCAP::Cli::valid_index_by_factor('index', $opts{'index'}, $max, 1);
    }
  }
  elsif(exists $opts{'index'}) {
    die "ERROR: -index cannot be defined without -process\n";
  }

  # now safe to apply defaults
  $opts{'threads'} = 1 unless(defined $opts{'threads'});
  $opts{'seqtype'} = 'WGS' unless(defined $opts{'seqtype'});

  my $tmpdir = File::Spec->catdir($opts{'outdir'}, 'tmpPindel');
  make_path($tmpdir) unless(-d $tmpdir);
  my $progress = File::Spec->catdir($tmpdir, 'progress');
  make_path($progress) unless(-d $progress);
  my $logs = File::Spec->catdir($tmpdir, 'logs');
  make_path($logs) unless(-d $logs);

  $opts{'tmp'} = $tmpdir;

  return \%opts;
}

__END__

=head1 pindel.pl

Reference implementation of Cancer Genome Project indel calling
pipeline.

=head1 SYNOPSIS

pindel.pl [options]

  Required parameters:
    -outdir    -o   Folder to output result to.
    -reference -r   Path to reference genome file *.fa[.gz]
    -tumour    -t   Tumour BAM file
    -normal    -n   Normal BAM file
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
    -skipgerm  -sg  Don't output events with more evidence in normal BAM.
    -cpus      -c   Number of cores to use. [1]
                     - recommend max 4 during 'input' process.
    -softfil   -sf  VCF filter rules to be indicated in INFO field as soft flags

  Targeted processing (further detail under OPTIONS):
    -process   -p   Only process this step then exit, optionally set -index
    -index     -i   Optionally restrict '-p' to single job

  Other:
    -help      -h   Brief help message.
    -man       -m   Full documentation.
    -version   -v   Version

  File list can be full file names or wildcard, e.g.
    pindel.pl -c 4 -r some/genome.fa[.gz] -o myout -t tumour.bam -n normal.bam

  Run with '-m' for possible input file types.

=head1 OPTIONS

=over 2

=item B<-process>

Available processes for this tool are:

  input
  pindel
  pin2vcf
  merge

=item B<-index>

Possible index ranges for processes above are:

  input   = 1..2
  pindel  = 1..<total_refs_less_exclude>
  pin2vcf = 1..<total_refs_less_exclude>
  merge   = 1
  flag    = 1

If you want STDOUT/ERR to screen ensure index is set even for single job steps.

=back
