#!/usr/bin/env perl

use strict;
use warnings FATAL => 'all';
use autodie qw(:all);
use Cwd qw(abs_path);
use File::Path qw(make_path);
use File::Spec::Functions;
use FindBin qw($Bin);
use Getopt::Long;
use lib "$Bin/../lib";
use Pod::Usage qw(pod2usage);

use PCAP::Cli;
use Sanger::CGP::Pindel::OutputGen::VcfBlatAugment;

{
  my $options = setup();
  my $augment = Sanger::CGP::Pindel::OutputGen::VcfBlatAugment->new(
    input => $options->{input},
    ref => $options->{ref},
    ofh => $options->{output},
    sam => $options->{align},
    hts_files => $options->{hts},
    outpath => $options->{outpath},
  );

  $augment->output_header;
  $augment->process_records;
}


sub setup{
  my %opts = (
    'cmd' => join(" ", $0, @ARGV),
    'hts' => [],
  );
  my @hts_files;
  GetOptions( 'h|help' => \$opts{h},
              'm|man' => \$opts{m},
              'v|version' => \$opts{v},
              'o|output=s' => \$opts{output},
              'r|ref=s' => \$opts{ref},
              'i|input=s' => \$opts{input},
              'd|debug' => \$opts{debug},
              'hts=s@' => \@hts_files,
  );

  $opts{hts} = [split(/,/,join(',',@hts_files))];


  if(defined $opts{'v'}) {
    printf "Version: %s\n", Sanger::CGP::Pindel::Implement->VERSION;
    exit;
  }

  pod2usage(-verbose => 1) if(defined $opts{h});
  pod2usage(-verbose => 2) if(defined $opts{m});

  PCAP::Cli::file_for_reading('ref', $opts{ref});
  PCAP::Cli::file_for_reading('input', $opts{input});
  for my $t(@{$opts{hts}}) {
    PCAP::Cli::file_for_reading('hts', $t);
    unless(-e $t.'.bai' || -e $t.'.csi' || -e $t.'.crai') {
      die "ERROR: Unable to find appropriate index file for $t\n";
    }
    unless(-e $t.'.bas') {
      die "ERROR: Unable to find *.bas file for $t\n";
    }
  }

  $opts{outpath} = $opts{output};
  make_path($opts{outpath}) unless(-e $opts{outpath});

  $opts{align} = catfile($opts{outpath}, 'data.sam');

  open my $ofh, '>', catfile($opts{outpath}, 'data.vcf');
  $opts{output} = $ofh;

  return \%opts;
}

__END__

=head1 NAME

pindel_blat_vaf.pl - Takes a raw Pindel VCF and bam file to add accurate counts.

=head1 SYNOPSIS

pindel_blat_vaf.pl [options]

  Required parameters:
    -ref       -r   File path to the reference file used to provide the coordinate system.
    -input     -i   VCF file to read in.
    -hts            BAM/CRAM file for associated sample.
    -output    -o   Directory for VCF output (gz compressed) and colocated sample bams

  Other:
    -help      -h   Brief help message.
    -man       -m   Full documentation.
    -version   -v   Prints the version number.

=head1 DESCRIPTION

B<pindel_blat_vaf.pl> will attempt to generate a vcf with expanded counts and VAF.

For every variant called by Pindel a blat will be performed and the results merged into a single vcf record.

=cut
