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
