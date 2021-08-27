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
use warnings FATAL=>'all';
use autodie;

use IO::Uncompress::Gunzip qw(gunzip $GunzipError);

if(@ARGV < 2) {
  warn "\nUSAGE: $0 normal_panel.gff3[.gz] sampleToRemove_1 sampleToRemove_2 [...]\n\n";
  warn "Output is to stdout, recommended usage:\n";
  warn "\t$0 normal_panel.gff3[.gz] sampleToRemove_1 sampleToRemove_2 [...] | bgzip -c > new_panel.gff3.gz\n";
  exit(1);
}

my $input = shift @ARGV;
my @samples = @ARGV;

my $z = IO::Uncompress::Gunzip->new($input, MultiStream => 1) or die "gunzip failed: $GunzipError\n";
while(my $line = <$z>) {
  if($. < 3) {
    print $line;
    next;
  }
  if($. == 3) {
    print clean_header($line, \@samples);
    next;
  }
  print clean_record($line, \@samples);
}
close $z;

sub clean_record {
  my ($line, $samples) = @_;
  chomp $line;
  my @columns = split /\t/, $line;
  my ($summary, @sample_stat) = split ',', pop @columns;
  my $found = 0;
  my @cleaned;

  for my $ss(@sample_stat) {
    my $matched = 0;
    for my $samp (@{$samples}) {
      if($ss =~ m/^$samp=/) {
        $matched++;
        next;
      }
    }
    if($matched > 0){
      $found++;
    }
    else {
      push @cleaned, $ss;
    }
  }
  if($found == 0) {
    return $line."\n";
  }
  my ($orig_count) = $summary =~ m/(\d+)$/;
  my $fixed_count = $orig_count - $found;
  return q{} if($fixed_count == 0); # removes line completely
  my $corrected = join "\t", @columns, join ',', 'SAMPLE_COUNT='.$fixed_count, @cleaned;
  return $corrected."\n";
}

sub clean_header {
  my ($line, $samples) = @_;
  chomp $line;
  my @elements = split /\s/, $line;
  my @new;
  OUTER: for my $e (@elements) {
    for my $samp (@{$samples}) {
      next OUTER if $e =~ m/$samp/;
    }
    push @new, $e;
  }
  my $final = join q{ }, @new;
  $final .= "\n";
  return $final;
}
