#!/usr/bin/env perl

use strict;
use warnings FATAL => 'all';
use autodie qw(:all);
use Cwd qw(abs_path);
use Pod::Usage qw(pod2usage);
use FindBin qw($Bin);
use lib "$Bin/../lib";
use Getopt::Long;

{
  my $options = setup();
}

sub setup {
  my %opts = (
    suffix => '.vaf',
  );
  GetOptions( 'h|help' => \$opts{h},
              'm|man' => \$opts{'m'},
              'v|version' => \$opts{v},
              'i|input=s' => \$opts{input},
              's|suffix:s' => \$opts{suffix},
  );

  if(defined $opts{v}) {
    printf "Version: %s\n", Sanger::CGP::Pindel::Implement->VERSION;
    exit;
  }
  pod2usage(-verbose => 1) if(defined $opts{h});
  pod2usage(-verbose => 2) if(defined $opts{'m'});

  PCAP::Cli::file_for_reading('input', $opts{input});

  my @htsfiles = @ARGV;
  for my $hts(@htsfiles) {
    PCAP::Cli::file_for_reading('bam/cram files', $hts);
  }
  $opts{'htsfiles'} = \@htsfiles;

  return \%opts;

}

__END__

=head1 NAME

pindelCohortVafFill.pl - Takes a VCF and adds VAF for sample/event with no call.

=head1 SYNOPSIS

pindelCohortVafFill.pl [options] SAMPLE_1.bam SAMPLE_2.bam [...]

  SAMPLE.bam should have co-located *.bai and *.bas files.

  Required parameters:
    -input     -i   VCF file to read in.

  Optional parameters:
    -suffix    -s   Suffix for updated file [.vaf]

  Other:
    -help      -h   Brief help message.
    -man       -m   Full documentation.
    -version   -v   Prints the version number.

=head1 DESCRIPTION

B<pindelCohortVafFill.pl> Fills in VAF for sample/event combinations where no call was made.

There must be a BAM/CRAM for every sample indicated by the VCF header.

=cut
