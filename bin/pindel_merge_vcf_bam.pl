#!/usr/bin/perl

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

#Sanger::CGP::Pindel::OutputGen::BamUtil::sam_to_sorted_bam();

{
  my $options = setup();
  Sanger::CGP::Pindel::OutputGen::BamUtil::sam_to_sorted_bam($options->{'out'}.'_mt', $options->{'mt'});
  Sanger::CGP::Pindel::OutputGen::BamUtil::sam_to_sorted_bam($options->{'out'}.'_wt', $options->{'wt'});
  merge_vcf($options->{'out'}, $options->{'vcf'})
}

sub merge_vcf {
  my ($path_prefix, $vcf_files) = @_;
  my $new_vcf = $path_prefix.'.vcf';
  system(qq{grep '^#' $vcf_files->[0] > $new_vcf});
  system(qq{cat @{$vcf_files} | grep -v '^#' | sort -k1,1 -k2,2n >> $new_vcf});

  my $vcf_gz = $new_vcf.'.gz';
  my $command = which('bgzip');
  $command .= sprintf ' -c %s > %s', $new_vcf, $vcf_gz;
  system($command);

  $command = which('tabix');
  $command .= sprintf ' -p vcf %s', $vcf_gz;
  system($command);
  return 1;
}

sub setup {
  my %opts;
  $opts{'cmd'} = join " ", $0, @ARGV;
  GetOptions( 'h|help' => \$opts{'h'},
              'm|man' => \$opts{'m'},
              'o|out=s' => \$opts{'out'},
  ) or pod2usage(2);

  my $version = Sanger::CGP::Pindel::OutputGen->VERSION;

  if(defined $opts{'v'}){
    print "Version: $version\n";
    exit;
  }

  pod2usage(-verbose => 1) if(defined $opts{'h'});
  pod2usage(-verbose => 2) if(defined $opts{'m'});

  # the files from the command line pattern match
  my @bad_files;
  for my $file(@ARGV) {
    if($file =~ m/\.vcf$/) {
      push @{$opts{'vcf'}}, $file;
    }
    elsif($file =~ m/_([mw]t)\.sam$/) {
      push @{$opts{$1}}, $file;
    }
    else {
      push @bad_files, $file;
    }
  }

  die "ERROR: The following unexpected files were presented: \n\t".(join "\n\t", @bad_files)."\n"
    if(scalar @bad_files);

  return \%opts;
}

__END__

=head1 NAME

pindel_merge_vcf_bam.pl - Merges provided VCF and SAM files into final result files.

=head1 SYNOPSIS

pindel_merge_vcf_bam.pl [options] files...

  Input files must follow the expected naming convention of:
    *.vcf    - the unmerged VCF files
    *_mt.sam - the remapped reads from the tumour sample
    *_wt.sam - the remapped reads from the normal sample

  This matches the standard output of pindel_2_combined_vcf.pl.

  Required parameters:
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

  Other:
    -help      -h   Brief help message.
    -man       -m   Full documentation.
    -version   -v   Prints the version number.

  Example:
   pindel_merge_vcf_bam.pl -o someloc/tum_vs_norm in/*.vcf in/*_mt.sam in/*_wt.sam

