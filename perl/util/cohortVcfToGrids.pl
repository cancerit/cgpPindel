#!/usr/bin/perl
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
