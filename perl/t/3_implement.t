########## LICENCE ##########
# Copyright (c) 2014-2022 Genome Research Ltd.
#
# Author: CASM/Cancer IT <cgphelp@sanger.ac.uk>
#
# This file is part of cgpPindel.
#
# cgpPindel is free software: you can redistribute it and/or modify it under
# the terms of the GNU Affero General Public License as published by the Free
# Software Foundation; either version 3 of the License, or (at your option) any
# later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
# details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
########## LICENCE ##########


use strict;
use Data::Dumper;
use Test::More tests => 1;

use Sanger::CGP::Pindel::Implement;

use FindBin qw($Bin);
my $lib_path = "$Bin/../lib";
my $test_data_path = "$Bin/../t/data/";

my $reference = $test_data_path.'genome_excludeseq.fa';
my $exclude_str = 'MT,GL%,hs37d5,NC_007605';
my $exclude_file_regex = $test_data_path.'test_exclude.txt';
my @exp_seqs = qw( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y);


subtest 'Test valid_seqs' => sub {
  my %options = ('reference' => $reference,
                 'exclude' => $exclude_str);

  my @good_seqs = Sanger::CGP::Pindel::Implement::valid_seqs(\%options);
  is_deeply(\@good_seqs,\@exp_seqs,'Exclude string');

  %options = ('reference' => $reference,
               'excludef' => $exclude_file_regex,);
  @good_seqs = Sanger::CGP::Pindel::Implement::valid_seqs(\%options);
  is_deeply(\@good_seqs,\@exp_seqs,'Exclude file');
  done_testing();
};
