#!/usr/bin/perl

use strict;
use warnings FATAL=>'all';
use autodie;
use Data::Dumper;

use IO::Uncompress::Gunzip qw(gunzip $GunzipError);

my $input = $ARGV[0];

my $z = IO::Uncompress::Gunzip->new($input, MultiStream => 1) or die "gunzip failed: $GunzipError\n";
while(my $line = <$z>) {
  if($line =~ m/^#/) {
    next if(index($line, 'pindel_np_from_vcf.pl') != -1);
    print $line;
    next;
  }
  chomp $line;
  next if($line eq q{});
  my @dataset = split /\t/, $line;
  my $sample_data = pop @dataset;
  my @sample_set = split /,/, $sample_data;
  my $new_samples = shift @sample_set;
  my $s_count = 1;
  for(@sample_set) {
    my ($raw_samp, $counts) = $_ =~ m/([^=]+)(.+)/;
    $new_samples .= ','.($s_count++).'_'.$counts;
  }
  print join("\t", @dataset, $new_samples),"\n";
}
close $z;
