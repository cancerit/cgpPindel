#!/usr/bin/perl

use strict;
use IO::Uncompress::Gunzip qw($GunzipError) ;
use List::Util qw(sum);

my ($combined_vcf, $exclude_sample, $output) = @ARGV;
# use a basic bed style output so we can intersect etc

my ($wtp, $wtn, $mtp, $mtn) = (5,6,8,9);

my $z = new IO::Uncompress::Gunzip $combined_vcf, MultiStream => 1 or die "IO::Uncompress::Gunzip failed: $GunzipError\n";
open my $O_DPTH, '>', $output.'.depth' or die "Failed to create $output.depth";
open my $O_ALT, '>', $output.'.alt' or die "Failed to create $output.alt";
my @sample_order;
while (my $l = <$z>) {
  next if $l =~ m/^##/;
  chomp $l;
  my @data = split "\t", $l;
  my ($chr, $pos, $ref, $alt, $gt) = @data[0,1,3,4,8];
  my @samples = @data[9..$#data];
  # header junk
  if($chr eq '#CHROM'){
    @sample_order = @samples;
    print $O_DPTH 'ID';
    print $O_ALT 'ID';
    for my $idx (0..$#samples) {
      if($sample_order[$idx] eq $exclude_sample) {
        next;
      }
      print $O_DPTH "\t".$sample_order[$idx];
      print $O_ALT "\t".$sample_order[$idx];
    }
    print $O_DPTH "\n";
    print $O_ALT "\n";

    next;
  }
  # body
  $chr =~ s/^chr//;
  my $id = sprintf '%s_%d_%s_%s', $chr, $pos, $ref, $alt;
  my @total_depth = ( $id );
  my @alt_depth = ( $id );
  for my $idx (0..$#samples) {
    if($sample_order[$idx] eq $exclude_sample) {
      next;
    }
    my @gt_data = split ':', $samples[$idx];
    push @total_depth, sum (@gt_data[$wtp, $wtn, $mtp, $mtn]);
    push @alt_depth, sum (@gt_data[$mtp, $mtn]);
    #print "$gt : $samples[$idx] : $total_depth : $mut_depth\n";
  }
  print $O_DPTH join("\t", @total_depth)."\n";
  print $O_ALT join("\t", @alt_depth)."\n";
}
close $O_DPTH;
close $O_ALT;
close $z;
