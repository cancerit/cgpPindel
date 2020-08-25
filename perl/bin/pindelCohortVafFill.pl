#!/usr/bin/env perl

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
  my $threads = PCAP::Threaded->new($options->{'threads'});

  if(!exists $options->{'process'} || $options->{'process'} eq 'split') {
    Sanger::CGP::Pindel::Implement::cohort_split($options);
  }

  if(!exists $options->{'process'} || $options->{'process'} eq 'fill') {
    for my $f(glob(catfile($options->{'split_dir'}, '*.vcf.gz'))) {
      push @{$options->{split_files}}, $f if($f =~ m{/\d+.vcf.gz$});
    }
    $threads->add_function('fill', \&Sanger::CGP::Pindel::Implement::fill_split_vaf);
    $threads->run(scalar @{$options->{split_files}}, 'fill', $options);
  }

  if(!exists $options->{'process'} || $options->{'process'} eq 'bams') {
    $threads->add_function('bams', \&Sanger::CGP::Pindel::Implement::merge_vaf_bams);
    $threads->run(scalar @{$options->{primary_hts}}, 'bams', $options);
  }

  if(!exists $options->{'process'} || $options->{'process'} eq 'finalise') {
    Sanger::CGP::Pindel::Implement::fill_vcf_merge($options);
    move(catdir($options->{tmp}, 'logs'), catdir($options->{output}, 'logs'));
    remove_tree($options->{tmp});
  }


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
              'c|cpus:i' => \$opts{'threads'},
              'p|process:s' => \$opts{'process'},
              'i|index:i' => \$opts{'index'},
              'l|limit:i' => \$opts{'limit'},
              'a|abort' => \$opts{'abort'},
              'd|data=s' => \$opts{'data'},
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

  for my $param(qw(input ref output)) {
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

  $opts{tmp} = catdir($opts{output}, 'tmpCohortVafFill');
  make_path($opts{tmp});

  $opts{split_dir} = catdir($opts{tmp}, 'split');
  make_path($opts{split_dir});
  $opts{fill_dir} = catdir($opts{tmp}, 'fill');
  make_path($opts{fill_dir});

  make_path(catdir($opts{tmp}, 'logs'));

  return \%opts;
}

__END__

=head1 NAME

pindelCohortVafFill.pl - Takes merged cohort VCF and fills in gaps in farm friendly manner.

=head1 SYNOPSIS

pindelCohortVafFill.pl [options] -i ... -o ... -r ... -d ...

  Required parameters:
    -file      -f   VCF file to read in.
    -output    -o   Workspace directory and final output.
    -ref       -r   File path to the reference file used to provide the coordinate system.
    -data      -d   File containing list of sequence data files for all samples used in "-file"
                    - format: tab separated BWA mapping followed by pindel_cohort reads, one sample per line.

                        sample_A_bwa.bam<TAB>sample_A_pindel.bam
                        sample_B_bwa.bam<TAB>sample_B_pindel.bam

  Optional parameters:
    -name      -n   Stub name for final output files [$output/cohort...]
    -size      -s   Number of Sample/event combinations per-file when processing [10000]
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

=head1 DESCRIPTION

B<pindelCohortVafFill.pl> Takes merged cohort VCF and fills in gaps in farm friendly manner.

One additional file is generated containg events where no fill-in is required.

Each original sample BAM/CRAM (normally BWA mapping) needs to be paired with the one generated by pindelCohort.pl

Each sample BAM/CRAM also needs colocated index and bas file.

=cut
