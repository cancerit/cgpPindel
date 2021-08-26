#!/usr/bin/perl
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
