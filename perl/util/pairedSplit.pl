#!/usr/bin/env perl

use strict;
use v5.16;
use Data::Dumper;

my ($t_name, $n_name, $in) = @ARGV;
open my $IF, '<', $in or die "Failed to open $in: $!";
my ($t_pos, $n_pos);
while(my $l = <$IF>) {
  next if ($l =~ m/^##/);
  chomp $l;
  my @bits = split /\t/, $l;
  if($l =~ m/^#CHRO/) {
    ($t_pos, $n_pos) = find_samples(\@bits, $t_name, $n_name);
    next;
  }
  process_record($bits[$t_pos], $bits[$n_pos], $l);
}
close $IF;

sub process_record {
  my ($t_geno, $n_geno, $record) = @_;
  if($n_geno eq q{.}) {
    return;
  }
  if($t_geno eq q{.}) {
    return;
  }
  else {
    my (undef, $s1, $s2, $pp, $np, $wp, $wn, $wm, $mp,$mn, $mm, $vaf) = split ':', $t_geno;
    return if(($mp + $mn + $wp + $wn) > 60); # example is 60x
    return if($mm >= 0.035);
    return if($vaf < 0.05);
    return if($mp == 0);
    return if($mn == 0);
    return if($wp == 0);
    return if($wn == 0);
    say $record;
    #if($mm < 0.035 && $vaf >= 0.02 && $mp > 0 && $mn > 0 && ) {
    #  say $record;
    #}
  }

}

sub find_samples {
  my ($bits, $t_name, $n_name) = @_;
  my ($t_pos, $n_pos);
  my $a_s = @{$bits} - 1;
  for my $i(0..$a_s) {
    if($bits->[$i] eq $t_name) {
      $t_pos = $i;
    }
    if($bits->[$i] eq $n_name) {
      $n_pos = $i;
    }
  }
  unless(defined $t_pos) {
    die "Failed to find $t_name\n"
  }
  unless(defined $n_pos) {
    die "Failed to find $n_name\n"
  }
  return ($t_pos, $n_pos);
}
