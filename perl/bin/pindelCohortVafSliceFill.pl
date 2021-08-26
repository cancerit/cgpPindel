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
use File::Path qw(make_path);
use File::Spec::Functions;
use FindBin qw($Bin);
use Getopt::Long;
use IO::Compress::Gzip qw(:constants gzip $GzipError);
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
    hts_files => $options->{hts_files},
    outpath => $options->{outpath},
    fill_in => 1,
  );

  $augment->output_header;
  $augment->process_records;
  $augment->sam_to_bam;
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
              'd|data=s' => \$opts{data},
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

  my @htsfiles;
  open my $D, '<', $opts{data};
  while(my $hts_file = <$D>) {
    chomp $hts_file;
    PCAP::Cli::file_for_reading('bam/cram files', $hts_file);
    push @htsfiles, $hts_file;
  }
  close $D;

  $opts{hts_files} = \@htsfiles;

  $opts{outpath} = $opts{output};
  make_path($opts{outpath}) unless(-e $opts{outpath});

  my $ofh = new IO::Compress::Gzip catfile($opts{output}, 'slice.vcf.gz'), -Level => Z_BEST_SPEED or die "IO::Compress::Gzip failed: $GzipError\n";
  $opts{output} = $ofh;

  return \%opts;
}

__END__

=head1 NAME

pindelCohortVafSliceFill.pl - Takes a VCF and adds VAF for sample/event with no call.

=head1 SYNOPSIS

pindelCohortVafSliceFill.pl [options]

  Required parameters:
    -ref       -r   File path to the reference file used to provide the coordinate system.
    -input     -i   VCF file to read in.
    -output    -o   Directory for VCF output (gz compressed) and colocated sample bams
    -data      -d   File containing list of BWA mappingfiles for all samples used in "-input"
                    - format: one BWA bam/cram file per line, expects co-located *.bai

  Other:
    -help      -h   Brief help message.
    -man       -m   Full documentation.
    -version   -v   Prints the version number.

=head1 DESCRIPTION

B<pindelCohortVafSliceFill.pl> Fills in VAF for sample/event combinations where no call was made.

There must be a BAM/CRAM for every sample indicated by the VCF header.

=cut
