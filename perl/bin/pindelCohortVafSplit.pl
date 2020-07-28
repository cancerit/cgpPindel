#!/usr/bin/env perl

use strict;
use warnings FATAL => 'all';
use autodie qw(:all);
use Cwd qw(abs_path);
use Pod::Usage qw(pod2usage);
use FindBin qw($Bin);
use lib "$Bin/../lib";
use Getopt::Long;
use File::Path qw(make_path);
use File::Spec::Functions;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use PCAP::Cli;

{
  my $options = setup();
  split_data($options->{input}, $options->{output}, $options->{size});
}

sub split_data {
  my ($input, $outdir, $max_e) = @_;
  make_path($outdir);
  my $part_fmt = '%s/%04d.vcf';
  my @header;
  my $no_vaf_count = 0;
  my $no_vaf_file_no = 0;

  my $complete_rec = catfile($outdir, 'complete_rec.vcf');
  my $no_vaf_file = sprintf $part_fmt, $outdir, $no_vaf_file_no++;

  open my $COMP, '>', $complete_rec;
  my $NO_VAF;

  my $z = IO::Uncompress::Gunzip->new($input, MultiStream => 1) or die "gunzip failed: $GunzipError\n";
  while(my $line = <$z>) {
    if($line =~ m/^#/) {
      push @header, $line;
      print $COMP $line;
      next;
    }
    unless($NO_VAF) {
      warn "Creating $no_vaf_file\n";
      open $NO_VAF, '>', $no_vaf_file;
      print $NO_VAF join q{}, @header;
    }
    chomp $line;
    # col 9+ containing '.' means work to be done
    my ($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @samples) = split /\t/, $line;
    my $n_vaf = 0;
    map { if($_ eq q{.}) {$n_vaf++} } @samples;
    if($n_vaf) {
      print $NO_VAF $line."\n";
      $no_vaf_count += $n_vaf;
      if($no_vaf_count >= $max_e) {
        close $NO_VAF;
        $no_vaf_file = sprintf $part_fmt, $outdir, $no_vaf_file_no++;
        warn "\nCreating $no_vaf_file\n";
        open $NO_VAF, '>', $no_vaf_file;
        print $NO_VAF join q{}, @header;
        $no_vaf_count = 0;
      }
    }
    else {
      print $COMP $line."\n";
    }
  }
  close $COMP;
  close $NO_VAF;
}

sub setup {
  my %opts = (
    'size' => 10000,
  );
  GetOptions( 'h|help' => \$opts{h},
              'm|man' => \$opts{'m'},
              'v|version' => \$opts{v},
              'i|input=s' => \$opts{input},
              'o|output=s' => \$opts{output},
              's|size:i' => \$opts{size},
  );

  if(defined $opts{v}) {
    printf "Version: %s\n", Sanger::CGP::Pindel::Implement->VERSION;
    exit;
  }
  pod2usage(-verbose => 1) if(defined $opts{h});
  pod2usage(-verbose => 2) if(defined $opts{'m'});

  PCAP::Cli::file_for_reading('input', $opts{input});

  return \%opts;
}

__END__

=head1 NAME

pindelCohortVafSplit.pl - Takes merged cohort VCF and splits into even sized files for VAF fill-in.

=head1 SYNOPSIS

pindelCohortVafSplit.pl [options]

  Required parameters:
    -input     -i   VCF file to read in.
    -output    -o   Directory for split VCFs

  Optional parameters:
    -size      -s   Number of Sample/event combinations per-file [10000]

  Other:
    -help      -h   Brief help message.
    -man       -m   Full documentation.
    -version   -v   Prints the version number.

=head1 DESCRIPTION

B<pindelCohortVafSplit.pl> will generate a set of split files for VAF fill-in processing.

One additional file is generated containg events where no fill-in is required.

=cut
