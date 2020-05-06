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
use File::Temp qw(tempdir);
use Test::More;
use Test::Fatal;
use File::Spec::Functions;
use Const::Fast qw(const);
use FindBin qw($Bin);

const my $MODULE => 'Sanger::CGP::Pindel::OutputGen::VcfBlatAugment';
const my $DATA => "$Bin/data/blat";
const my @HEADER_ENDS => do {
    no warnings 'qw';
    qw(#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT PD26988a);
};
const my $HEADER_LINES => 3391;

const my $DATA_ARR_D => [qw(chr10 11201 id CA C 390 . PC=D;RS=11201;RE=11205;LEN=1;S1=11;S2=849.236;REP=3 GT:PP:NP ./.:10:0)];
const my $RES_ARR_D => [13, 3, 9, 3, 0.429];
const my $DATA_ARR_DI => [qw(chr10 22777 id AGAAACTGTG ACTGTGAGATAGATATATATAGATAGATATAT 105 . PC=DI;RS=22777;RE=22787;LEN=9;S1=6;REP=0 GT:PP:NP ./.:0:5)];
const my $RES_ARR_DI => [2, 3, 2, 9, 0.688];
const my $DATA_ARR_SI => [qw(chr10 11643 id C CG 150 . PC=I;RS=11643;RE=11649;LEN=1;S1=6;S2=421.908;REP=4 GT:PP:NP ./.:0:5)];
const my $RES_ARR_SI => [7, 14, 7, 10, 0.447];

my ($stdout_fh, $buffer);

subtest 'Initialisation checks' => sub {
  use_ok($MODULE);
  new_vba(catfile($DATA, 'D.vcf'));
};

subtest 'Header checks' => sub {
  my $vba = new_vba(catfile($DATA, 'D.vcf'));
  ok($vba->output_header);
  my @lines = split /\n/, $buffer;
  is(scalar @lines, $HEADER_LINES, 'Expected number of header lines');
  is($lines[-1], join("\t", @HEADER_ENDS), 'Expected final header line');
};

subtest 'Simple Deletion checks' => sub {
  my $vba = new_vba(catfile($DATA, 'D.vcf'));
  my $v_h = $vba->to_data_hash($DATA_ARR_D);
  my @gt_out = $vba->blat_record($v_h);
  is_deeply(\@gt_out, $RES_ARR_D);
};

subtest 'Simple Insertion checks' => sub {
  my $vba = new_vba(catfile($DATA, 'SI.vcf'));
  my $v_h = $vba->to_data_hash($DATA_ARR_SI);
  my @gt_out = $vba->blat_record($v_h);
  is_deeply(\@gt_out, $RES_ARR_SI);
};

subtest 'Complex event checks' => sub {
  my $vba = new_vba(catfile($DATA, 'DI.vcf'));
  my $v_h = $vba->to_data_hash($DATA_ARR_DI);
  my @gt_out = $vba->blat_record($v_h);
  is_deeply(\@gt_out, $RES_ARR_DI);
};

done_testing();


sub new_vba {
  my $vcf = shift;
  return new_ok($MODULE, [
    input => $vcf,
    ref => catfile($DATA, 'chr10_1-23700.fa'),
    ofh => buffer_fh(),
    hts_file => catfile($DATA, 'test.bam'),
  ]);
}

sub buffer_fh {
  if(defined $stdout_fh) {
    close $stdout_fh;
  }
  $buffer = q{};
  open $stdout_fh, ">", \$buffer or die $!;
  return $stdout_fh;
}



__END__
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

