#!/usr/bin/perl

########## LICENCE ##########
# Copyright (c) 2014-2021 Genome Research Ltd.
#
# Author: CASM/Cancer IT <cgphelp@sanger.ac.uk>
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

use File::Which qw(which);
use Getopt::Long;
use Try::Tiny;
use Carp;
use Pod::Usage qw(pod2usage);

use Sanger::CGP::Pindel::OutputGen::BamUtil;
use Sanger::CGP::Pindel::Implement;

{
  my $options = setup();
  Sanger::CGP::Pindel::OutputGen::BamUtil::sam_to_sorted_bam($options->{'out'}.'_mt',
                                                             $options->{'indir'},
                                                             $options->{'mt'},
                                                             $options->{'cram'},
                                                             $options->{'csi'},
                                                             $options->{'ref'},
                                                             );
  Sanger::CGP::Pindel::OutputGen::BamUtil::sam_to_sorted_bam($options->{'out'}.'_wt',
                                                             $options->{'indir'},
                                                             $options->{'wt'},
                                                             $options->{'cram'},
                                                             $options->{'csi'},
                                                             $options->{'ref'},
                                                             );
  merge_vcf($options->{'out'}, $options->{'indir'}, $options->{'vcf'})
}

sub merge_vcf {
  my ($path_prefix, $indir, $vcf_files) = @_;
  $vcf_files = Sanger::CGP::Pindel::Implement::fragmented_files($indir, $vcf_files, '#', 'FINAL_MERGED.vcf');
  my $new_vcf = $path_prefix.'.vcf';
  system(qq{cd $indir; grep -B 100000000 -m 1 '^#CHROM' $vcf_files->[0] > $new_vcf});
  system(qq{cd $indir; cat @{$vcf_files} | grep -v '^#' | sort -k1,1 -k2,2n -k4,4 -k 5,5 >> $new_vcf});

  my $vcf_gz = $new_vcf.'.gz';
  my $command = which('bgzip');
  $command .= sprintf ' -c %s > %s', $new_vcf, $vcf_gz;
  system($command);

  $command = which('tabix');
  $command .= sprintf ' -p vcf %s', $vcf_gz;
  system($command);
  unlink $new_vcf; # remove the VCF you created prior to compression
  return 1;
}

sub setup {
  my %opts;
  $opts{'cmd'} = join " ", $0, @ARGV;

  GetOptions( 'h|help' => \$opts{'h'},
              'm|man' => \$opts{'m'},
              'v|version' => \$opts{'v'},
              'o|out=s' => \$opts{'out'},
              'i|indir=s' => \$opts{'indir'},
              'c|cram' => \$opts{'cram'},
              'r|ref:s' => \$opts{'ref'},
              's|csi' => \$opts{'csi'},
  ) or pod2usage(2);

  if(defined $opts{'v'}){
    print 'Version: ',Sanger::CGP::Pindel::OutputGen::BamUtil->VERSION,"\n";
    exit;
  }

  pod2usage(-verbose => 1) if(defined $opts{'h'});
  pod2usage(-verbose => 2) if(defined $opts{'m'});

  # the files from the command line pattern match
#  my @bad_files;
#  for my $file(@ARGV) {

  $opts{'indir'} .= '/' if($opts{'indir'} !~ m|/$|);
  opendir(my $dh, $opts{'indir'});
  while(my $file = readdir $dh) {
    # don't store the full_path, as it makes the memory footprint too large
    my $full_path = $opts{'indir'}.$file;
    next unless(-f $full_path);
    if($file =~ m/\.vcf$/) {
      next if($file eq 'FINAL_MERGED.vcf');
      push @{$opts{'vcf'}}, $file;
    }
    elsif($file =~ m/_([mw]t)\.sam(.gz)?$/) {
      next if($file =~ m/^FINAL_MERGED\.sam.gz$/); # this is transient for each type
      push @{$opts{$1}}, $file;
    }
  }
  closedir($dh);

#  die "ERROR: The following unexpected files were presented: \n\t".(join "\n\t", @bad_files)."\n"
#    if(scalar @bad_files);

warn scalar(@{$opts{'vcf'}})." vcf files\n";
warn scalar(@{$opts{'mt'}})." mt sam files\n";
warn scalar(@{$opts{'wt'}})." wt sam files\n";

  return \%opts;
}

__END__

=head1 NAME

pindel_merge_vcf_bam.pl - Merges provided VCF and SAM files into final result files.

=head1 SYNOPSIS

pindel_merge_vcf_bam.pl [options] files...

  Files found in '-indir' must follow the expected naming convention of:
    *.vcf    - the unmerged VCF files
    *_mt.sam - the remapped reads from the tumour sample
    *_wt.sam - the remapped reads from the normal sample

  This matches the standard output of pindel_2_combined_vcf.pl.

  Required parameters:
    -indir     -i   Directory containing per-contig files for:
                     - *.vcf
                     - *_mt.sam
                     - *_wt.sam
    -out       -o   Output stub for final files e.g.
                    somepath/sample_vs_sample, results in:

                      VCF+index
                        somepath/sample_vs_sample.vcf.gz
                        somepath/sample_vs_sample.vcf.gz.tbi

                      BAM+index for tumour/mutant reads
                        somepath/sample_vs_sample_mt.bam
                        somepath/sample_vs_sample_mt.bam.bai

                      BAM+index for normal/wildtype reads
                        somepath/sample_vs_sample_wt.bam
                        somepath/sample_vs_sample_wt.bam.bai

  Optional parameters:
    -csi        -s  Specify for csi index when BAM files generated
    -cram       -c  Specify for CRAM as final output (crai index)
    -ref        -r  genome.fa file, required for CRAM output

  Other:
    -help      -h   Brief help message.
    -man       -m   Full documentation.
    -version   -v   Prints the version number.

  Example:
   pindel_merge_vcf_bam.pl -o someloc/tum_vs_norm in/*.vcf in/*_mt.sam in/*_wt.sam
