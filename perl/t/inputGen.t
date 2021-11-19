########## LICENCE ##########
# Copyright (c) 2014-2021 Genome Research Ltd.
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
use File::Temp qw(tempdir);
use Test::More;
use Test::Fatal;
use Const::Fast qw(const);
use FindBin qw($Bin);

const my $MODULE => 'Sanger::CGP::Pindel::InputGen';
const my $DATA => "$Bin/data";

const my $RECORD_SET => [
  [ "\@2:2428:29677:5760/1_RG461144\nGTTAGGGTTAGGGTTAGGGTTGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAAGATAGGAAGAGCACA\n+\t22\t10001\t38\t364\tSAMPLE",
    "\@2:2428:29677:5760/2_RG461144\nTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCAACCCAAACCCAAACCCAAACA\n-\t22\t10130\t38\t364\tSAMPLE"
  ]
];
# length of each record plus line feed each (added in function)
const my $RECORD_OUT_BYTES => (length $RECORD_SET->[0]->[0]) + (length $RECORD_SET->[0]->[1]) +2;

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

subtest 'reads_to_disk checks' => sub{
  $obj = new_ok($MODULE,
                [ "$DATA/test.bam",
                  undef,
                  "$DATA/genome_22.fa"]);
  my $out_folder = tempdir( 'pindelTests_XXXX', CLEANUP => 1 );

  $obj->set_outdir($out_folder);
  ok($obj->reads_to_disk($RECORD_SET), 'create output');
  is($obj->{'rname_bytes'}->{'22'}, $RECORD_OUT_BYTES, 'Verify bytes written captured');

  ok($obj->validate, 'validate returns true')
};

done_testing();
