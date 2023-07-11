#!/usr/bin/env perl
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
use warnings FATAL => 'all';
use autodie qw(:all);
use File::Basename;
use Capture::Tiny qw(capture);
use Const::Fast qw(const);
use File::Temp;
use File::Path qw(make_path);
use File::Spec::Functions;

const my $USAGE => sprintf "USAGE: %s in.vcf <lines-per-split> <outputDir>\n", basename($0);
const my $SRT_COUNT => q{bash -c "set -o pipefail; (grep -B 100000  -m 1 '^#CHRO' %s &&  grep -v '^#' %s | sort -s -S 1G -k 1,1 -k 2,2n -k 4,4 -k 5,5) | tee %s | grep -cv '^#'"};
const my $NO_HEAD_SPLIT => q{bash -c "set -o pipefail; grep -v '^#' %s | split -a 4 --additional-suffix=.vcf -l %d - %s"};
const my $CAPTURE_HEADER => q{grep -B 100000  -m 1 '^#CHRO' %s};

if(@ARGV < 3) {
  die $USAGE;
}

my ($in_vcf, $split_lines, $outdir) = @ARGV;

die "ERROR: Absent or Empty file: $in_vcf" unless(-e $in_vcf && -s _ > 0);

make_path($outdir) unless(-e $outdir);

my $tmp_dir = File::Temp->newdir(DIR=> $outdir, CLEANUP => 1);
my $srt_vcf = catfile($tmp_dir, 'srt.vcf');

# first we sort and capture the number of variants
my $c_srt_count = sprintf $SRT_COUNT, $in_vcf, $in_vcf, $srt_vcf;
my ($c_out, $c_err, $c_exit) = capture { system($c_srt_count); };
if($c_exit > 1) { # allow 1 as could be 0 events to work with
  warn "An error occurred while executing $c_srt_count\n";
  warn "\tERROR$c_err\n";
  exit $c_exit;
}
chomp $c_out;
my $events = $c_out;
die "ERROR:  Did not get a count of events" if($events !~ m/^\d+$/);

my $prefix = catfile($tmp_dir, 'split_');
my $c_split = sprintf $NO_HEAD_SPLIT, $srt_vcf, $split_lines, $prefix;
($c_out, $c_err, $c_exit) = capture { system($c_split); };
if($c_exit > 1) { # allow 1 as could be 0 events to work with
  warn "An error occurred while executing $c_split\n";
  warn "\tERROR$c_err\n";
  exit $c_exit;
}

my $c_header = sprintf $CAPTURE_HEADER, $in_vcf;
($c_out, $c_err, $c_exit) = capture { system($c_header); };
if($c_exit > 0) {
  warn "An error occurred while executing $c_split\n";
  warn "\tERROR$c_err\n";
  exit $c_exit;
}
my $vcf_head = $c_out;

# now need to convert all the files to valid VCF:
opendir(my $dh, $tmp_dir) || die "Can't opendir $tmp_dir: $!";
while (readdir $dh) {
  next unless($_ =~ m/^split_[a-z]{4}\.vcf$/);
  my $split_vcf = catfile($tmp_dir, $_);
  my $vcf_out = catfile($outdir, $_);
  open my $ofh, '>', $vcf_out;
  print $ofh $vcf_head;
  close $ofh;
  system("cat $split_vcf >> $vcf_out");
}
closedir $dh;
