#!/usr/bin/perl

use strict;

my ($vaf_correct, $exclude_sample, $output) = @ARGV;

open my $O_VC, '<', $vaf_correct;

my $header = <$O_VC>;
chomp $header;
my @head_items = split "\t", $header;
my ($chr, $pos, $ref, $alt) = @head_items[2,3,4,5];
my %mtr_by_sample;
for my $i(6..$#head_items) {
  if($head_items[$i] =~ m/_MTR$/) {
    my $sample = $head_items[$i];
    $sample =~ s/_MTR$//;
    next if($sample eq $exclude_sample);
    $mtr_by_sample{$sample} = $i;
  }
}

my @samples = sort keys %mtr_by_sample;
my (@mtr_cols, @wtr_cols);
for my $s(@samples) {
  push @mtr_cols, $mtr_by_sample{$s};
  push @wtr_cols, $mtr_by_sample{$s} + 1;
}


open my $O_DPTH, '>', $output.'.depth' or die "Failed to create $output.depth";
open my $O_ALT, '>', $output.'.alt' or die "Failed to create $output.alt";

print $O_DPTH join("\t", "ID", @samples)."\n";
print $O_ALT join("\t", "ID", @samples)."\n";

while(my $l = <$O_VC>) {
  chomp $l;
  my @row = split "\t", $l;
  my ($chr, $pos, $ref, $alt) = @row[2,3,4,5];
  $chr =~ s/^chr//;
  my $id = sprintf '%s_%d_%s_%s', $chr, $pos, $ref, $alt;
  my @total_depth;
  my @alt_depth;
  for my $i(0..$#mtr_cols) {
    # same length (or broken)
    push @total_depth, $row[$mtr_cols[$i]] + $row[$wtr_cols[$i]];
    push @alt_depth, $row[$mtr_cols[$i]];
  }
  print $O_DPTH join("\t", $id, @total_depth)."\n";
  print $O_ALT join("\t", $id, @alt_depth)."\n";
}

close $O_DPTH;
close $O_ALT;
close $O_VC;
