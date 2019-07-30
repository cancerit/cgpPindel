########## LICENCE ##########
# Copyright (c) 2019 Genome Research Ltd.
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
use Test::More;
use Test::Fatal;
use Const::Fast qw(const);
use FindBin qw($Bin);

const my $MODULE => 'Sanger::CGP::Pindel::InputGen';
const my $DATA => "$Bin/data";

my $obj;
subtest 'Initialisation checks' => sub {
  use_ok($MODULE);
};

subtest 'corrupt_pindel_input checks' => sub {
  # File of correct size (no NUL)
  is( Sanger::CGP::Pindel::InputGen::corrupt_pindel_input("$DATA/inputGen-goodfile.txt.gz", 22),
      undef,
      'corrupt_pindel_input - well formed compressed file, expected size');
  # File of incorrect size (no NUL)
  is( Sanger::CGP::Pindel::InputGen::corrupt_pindel_input("$DATA/inputGen-goodfile.txt.gz", 9),
      "$DATA/inputGen-goodfile.txt.gz",
      'corrupt_pindel_input - well formed compressed file, UNexpected size');
  # File of incorrect size + NUL
  is( Sanger::CGP::Pindel::InputGen::corrupt_pindel_input("$DATA/inputGen-NUL.txt.gz", 22),
      "$DATA/inputGen-NUL.txt.gz",
      'corrupt_pindel_input - NUL character in compressed file, expected size');
  # File of correct size + NUL
  is( Sanger::CGP::Pindel::InputGen::corrupt_pindel_input("$DATA/inputGen-NUL.txt.gz", 9),
      "$DATA/inputGen-NUL.txt.gz",
      'corrupt_pindel_input - NUL character in compressed file, UNexpected size');
};

done_testing();

