#!/usr/bin/env perl
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
    debug => $options->{debug},
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
    -debug     -d   Turn on additional outputs.
    -help      -h   Brief help message.
    -man       -m   Full documentation.
    -version   -v   Prints the version number.

=head1 DESCRIPTION

B<pindel_blat_vaf.pl> will attempt to generate a vcf with expanded counts and VAF.

For every variant called by Pindel a blat will be performed and the results merged into a single vcf record.

=cut
