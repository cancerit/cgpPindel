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
const my @RM_PROCESS => qw(  process index limit exclude excludef badloci apid);
my %index_max = ( 'input'   => 2,
                  'pindel'  => -1,
                  'pin2vcf' => -1,
                  'merge'   => 1,
                  'flag'    => 1,);

{
  my $options = setup();
  Sanger::CGP::Pindel::Implement::prepare($options);
  my $threads = PCAP::Threaded->new($options->{'threads'});

  # register any process that can run in parallel here
  $threads->add_function('input', \&Sanger::CGP::Pindel::Implement::input, exists $options->{'index'} ? 1 : 2);
  $threads->add_function('pindel', \&Sanger::CGP::Pindel::Implement::pindel);
  $threads->add_function('pin2vcf', \&Sanger::CGP::Pindel::Implement::pindel_to_vcf);

  # start processes here (in correct order obviously), add conditions for skipping based on 'process' option
  $threads->run(2, 'input', $options) if(!exists $options->{'process'} || $options->{'process'} eq 'input');

  # count the valid input files, gives constant job count for downstream
  if(!exists $options->{'process'} || first { $options->{'process'} eq $_ } ('pindel', 'pin2vcf')) {
    my $jobs = Sanger::CGP::Pindel::Implement::determine_jobs($options); # method still needed to populate info
    $jobs = $options->{'limit'} if(exists $options->{'limit'} && defined $options->{'limit'});
    $threads->run($jobs, 'pindel', $options) if(!exists $options->{'process'} || $options->{'process'} eq 'pindel');
    $threads->run($jobs, 'pin2vcf', $options) if(!exists $options->{'process'} || $options->{'process'} eq 'pin2vcf');
  }

  Sanger::CGP::Pindel::Implement::merge_and_bam($options) if(!exists $options->{'process'} || $options->{'process'} eq 'merge');

  if(!exists $options->{'process'} || $options->{'process'} eq 'flag') {
    Sanger::CGP::Pindel::Implement::flag($options) unless($options->{'noflag'});
    cleanup($options) unless($options->{'debug'});
  }
}

sub cleanup {
  my $options = shift;
  my $tmpdir = $options->{'tmp'};
  move(File::Spec->catdir($tmpdir, 'logs'), File::Spec->catdir($options->{'outdir'}, 'logs')) || die $!;
  remove_tree $tmpdir if(-e $tmpdir);
  opendir(my $dh, $options->{'outdir'});
  unless($options->{'noflag'}) {
    while(readdir $dh) {
      unlink File::Spec->catfile($options->{'outdir'}, $_) if($_ =~ /\.vcf\.gz(\.tbi)?$/ && $_ !~ /\.flagged\.vcf\.gz(\.tbi)?$/);
    }
  }
  closedir $dh;
	return 0;
}

sub setup {
  my %opts;
  pod2usage(-msg  => "\nERROR: Option must be defined.\n", -verbose => 1,  -output => \*STDERR) if(scalar @ARGV == 0);
  $opts{'cmd'} = join " ", $0, @ARGV;
  GetOptions(
              'r|reference=s' => \$opts{'reference'},
              'o|outdir=s' => \$opts{'outdir'},
              't|tumour=s' => \$opts{'tumour'},
              'n|normal=s' => \$opts{'normal'},
              'e|exclude=s' => \$opts{'exclude'},
              'ef|exclude-file=s' => \$opts{'excludef'},
              'b|badloci=s' => \$opts{'badloci'},
              # these are specifically for pin2vcf
              'sp|species=s{0,}' => \@{$opts{'species'}},
              'as|assembly=s' => \$opts{'assembly'},
              'st|seqtype=s' => \$opts{'seqtype'},
              'sg|skipgerm' => \$opts{'skipgerm'},
              # specifically for FlagVCF
              's|simrep=s' => \$opts{'simrep'},
              'f|filters=s' => \$opts{'filters'},
              'g|genes=s' => \$opts{'genes'},
              'u|unmatched=s' => \$opts{'unmatched'},
              'sf|softfil=s' => \$opts{'softfil'},
              'a|apid:s' => \$opts{'apid'},
              # process management
              'c|cpus=i' => \$opts{'threads'},
              'l|limit=i' => \$opts{'limit'},
              'd|debug' => \$opts{'debug'},
              'p|process=s' => \$opts{'process'},
              'i|index=i' => \$opts{'index'},
              'noflag' => \$opts{'noflag'},
              # other
              'h|help' => \$opts{'h'},
              'm|man' => \$opts{'m'},
              'v|version' => \$opts{'version'},

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
  PCAP::Cli::file_for_reading('exclude file', $opts{'excludef'}) if(defined $opts{'excludef'});
  unless($opts{'noflag'}) {
    PCAP::Cli::file_for_reading('simrep', $opts{'simrep'});
    PCAP::Cli::file_for_reading('filters', $opts{'filters'});
    PCAP::Cli::file_for_reading('genes', $opts{'genes'});
    PCAP::Cli::file_for_reading('unmatched', $opts{'unmatched'});
    PCAP::Cli::file_for_reading('softfil', $opts{'softfil'}) if(defined $opts{'softfil'});
  }
  PCAP::Cli::out_dir_check('outdir', $opts{'outdir'});
  my $final_logs = File::Spec->catdir($opts{'outdir'}, 'logs');
  if(-e $final_logs) {
    warn "NOTE: Presence of '$final_logs' directory suggests successful complete analysis, please delete to rerun\n";
    exit 0;
  }

  for my $title (@RM_PROCESS){
    delete $opts{$title} unless(defined $opts{$title});
  }

  if(exists $opts{'process'}) {
    PCAP::Cli::valid_process('process', $opts{'process'}, \@VALID_PROCESS);
    if(exists $opts{'index'}) {
      my @valid_seqs = Sanger::CGP::Pindel::Implement::valid_seqs(\%opts);
      my $refs = scalar @valid_seqs;

      my $max = $index_max{$opts{'process'}};
      if($max==-1){
        if(exists $opts{'limit'}) {
          $max = $opts{'limit'} > $refs ? $refs : $opts{'limit'};
        }
        else {
      	  $max = $refs;
      	}
      }

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


  # make all things that appear to be paths complete (absolute not great if BAM/BAI in different locations)
  for my $key (keys %opts) {
    next unless( first {$key eq $_} qw(reference outdir tumour normal badloci simrep filters genes unmatched softfil));
    $opts{$key} = cwd().'/'.$opts{$key} if(defined $opts{$key} && -e $opts{$key} && $opts{$key} !~ m/^\//);
  }

  my $tmpdir = File::Spec->catdir($opts{'outdir'}, 'tmpPindel');
  make_path($tmpdir) unless(-d $tmpdir);
  my $progress = File::Spec->catdir($tmpdir, 'progress');
  make_path($progress) unless(-d $progress);
  my $logs = File::Spec->catdir($tmpdir, 'logs');
  make_path($logs) unless(-d $logs);

  $opts{'tmp'} = $tmpdir;

  if(scalar @{$opts{'species'}} > 0 ){
    $opts{'species'}="@{$opts{'species'}}";
  }
  else {
    delete $opts{'species'};
  }

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
    -tumour    -t   Tumour BAM/CRAM file (co-located index and bas files)
    -normal    -n   Normal BAM/CRAM file (co-located index and bas files)


  Flagging parameters:
    -noflag         Skip flagging
    -simrep    -s   Full path to tabix indexed simple/satellite repeats.
    -filter    -f   VCF filter rules file (see FlagVcf.pl for details)
    -genes     -g   Full path to tabix indexed coding gene footprints.
    -unmatched -u   Full path to tabix indexed gff3 or bed of unmatched normal panel
                      - gff3 = F010 & FF010
                      -  bed = FF021
                      - see also pindel_np_from_vcf.pl
    -softfil   -sf  Optional: Filter rules to be indicated in INFO field as soft flags
    -apid      -a   Optional: Analysis process ID (numeric)
                     - not necessary for external use

  Optional
    -seqtype        -st  Sequencing protocol, expect all input to match [WGS]
    -assembly       -as  Name of assembly in use
                           -  when not available in BAM header SQ line.
    -species        -sp  Species
                           -  when not available in BAM header SQ line.
    -exclude         -e  Exclude this list of ref sequences from processing, wildcard '%'
                           - comma separated, e.g. NC_007605,hs37d5,GL%
    -exclude-file    -ef File containing a list of ref sequences to exclude from processing, wildcard '%'
                           - can be used in place of the -e|-exclude option
    -badloci   -b   Tabix indexed BED file of locations to not accept as anchors
                     - e.g. hi-seq depth from UCSC
    -skipgerm  -sg  Don't output events with more evidence in normal BAM.
    -cpus      -c   Number of cores to use. [1]
                     - recommend max 4 during 'input' process.
    -limit     -l   When defined with '-cpus' internally thread concurrent processes.
                     - requires '-p', specifically for pindel/pin2vcf steps
    -debug     -d   Don't cleanup workarea on completion.


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
