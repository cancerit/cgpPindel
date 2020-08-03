#!/usr/bin/env perl

use strict;
use warnings FATAL => 'all';
use autodie qw(:all);
use Cwd qw(abs_path);
use Pod::Usage qw(pod2usage);
use FindBin qw($Bin);
use lib "$Bin/../lib";
use Getopt::Long;

use PCAP::Cli;
use Sanger::CGP::Pindel::OutputGen::VcfBlatAugment;

{
  my $options = setup();
  my $augment = Sanger::CGP::Pindel::OutputGen::VcfBlatAugment->new(
    input => $options->{input},
    ref => $options->{ref},
    ofh => $options->{output},
    sam => $options->{align},
    hts_files => $options->{hts_files},
    outpath => $options->{outpath},
    fill_in => 1,
  );

  $augment->output_header;
  $augment->process_records;

}

sub setup {
  my %opts = (
  );
  GetOptions( 'h|help' => \$opts{h},
              'm|man' => \$opts{'m'},
              'v|version' => \$opts{v},
              'i|input=s' => \$opts{input},
              'r|ref=s' => \$opts{ref},
              'o|output=s' => \$opts{output},
  );

  if(defined $opts{v}) {
    printf "Version: %s\n", Sanger::CGP::Pindel::Implement->VERSION;
    exit;
  }
  pod2usage(-verbose => 1) if(defined $opts{h});
  pod2usage(-verbose => 2) if(defined $opts{'m'});

  PCAP::Cli::file_for_reading('input', $opts{input});
  PCAP::Cli::file_for_reading('ref', $opts{ref});

  $opts{align} = $opts{output}.'.fill.sam' unless(defined $opts{align});

  my @htsfiles = @ARGV;
  for my $hts(@htsfiles) {
    PCAP::Cli::file_for_reading('bam/cram files', $hts);
  }
  $opts{'hts_files'} = \@htsfiles;

  $opts{outpath} = $opts{output};
  open my $ofh, '>', $opts{output};
  $opts{output} = $ofh;

  return \%opts;

}

__END__

=head1 NAME

pindelCohortVafSliceFill.pl - Takes a VCF and adds VAF for sample/event with no call.

=head1 SYNOPSIS

pindelCohortVafSliceFill.pl [options] SAMPLE_1.bam SAMPLE_2.bam [...]

  SAMPLE*.bam should have co-located *.bai and *.bas files.

  Required parameters:
    -ref       -r   File path to the reference file used to provide the coordinate system.
    -input     -i   VCF file to read in.
    -output    -o   File path for VCF output (not compressed), colcated sample bams

  Other:
    -help      -h   Brief help message.
    -man       -m   Full documentation.
    -version   -v   Prints the version number.

=head1 DESCRIPTION

B<pindelCohortVafSliceFill.pl> Fills in VAF for sample/event combinations where no call was made.

There must be a BAM/CRAM for every sample indicated by the VCF header.

=cut
