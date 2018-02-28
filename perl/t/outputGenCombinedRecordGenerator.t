########## LICENCE ##########
# Copyright (c) 2014-2018 Genome Research Ltd.
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
use FindBin qw($Bin);
use Data::Dumper;
use File::Spec;

use Bio::DB::HTS;

my $MODULE = 'Sanger::CGP::Pindel::OutputGen::CombinedRecordGenerator';

use_ok 'Sanger::CGP::Vcf::Sample';
use_ok 'Sanger::CGP::Vcf::Contig';
use_ok 'Sanger::CGP::Pindel::OutputGen::CombinedRecord';

my $test_data = "$Bin/data";
my $ref_file = File::Spec->catfile($test_data,'genome_22.fa');
my $pin_file = File::Spec->catfile($test_data,'test.txt_22');

my $mt_file = File::Spec->catfile($test_data,'test.bam');
my $wt_file = File::Spec->catfile($test_data,'test-BL.bam');

subtest 'Initialisation checks' => sub {
  use_ok($MODULE);

  my $fai = Bio::DB::HTS::Fai->load($ref_file);
  my $wt_sam = Bio::DB::HTS->new(-bam => $wt_file, -fasta => $ref_file);
  my $mt_sam = Bio::DB::HTS->new(-bam => $mt_file, -fasta => $ref_file);
  my $obj = new_ok($MODULE => [-path => $pin_file, -fai => $fai, -wt_sam => $wt_sam, -mt_sam => $mt_sam, -mutant_sample_name => 'test']); ## this will fail without a valid sam pindel file

  is($obj->{_wt_sam},$wt_sam,'Correctly setting wt sam object');
  is($obj->{_mt_sam},$mt_sam,'Correctly setting mt sam object');
  is($obj->{_mutant_sample_name},'test','Correctly setting sample name');

  eval{
    new $MODULE(-path => $pin_file, -fai => $fai, -wt_sam => $wt_sam, -mutant_sample_name => 'test');
  };unless($@){
  	fail('Initialisation checks: expected an exception to be thrown if no mt_sam object set.');
  }else{
  	pass('new method throws exception as expected if no mt_sam object arguments');
  }

  eval{
    new $MODULE(-path => $pin_file, -fai => $fai, -mt_sam => $mt_sam, -mutant_sample_name => 'test');
  };unless($@){
  	fail('Initialisation checks: expected an exception to be thrown if no wt_sam object set.');
  }else{
  	pass('new method throws exception as expected if no wt_sam object arguments');
  }

};



subtest 'Non-object funcions' => sub {


  subtest '_read_depth' => sub {

  	my $wt_sam = Bio::DB::HTS->new(-bam => $wt_file, -fasta => $ref_file);

  	my $expect = [{
  	  'HWUSI-EAS493_8289_FC30GNW_PE:6:86:991:775' => '1',
  	  'HWI-EAS255_8305_FC30G7J_PE:5:12:939:716' => '1',
  	  'HWI-EAS258_8310_FC30G9N_PE:5:81:148:430' => '1',
  	  'HWI-EAS107_8284_FC30GCE_PE:4:25:1041:38' => '1',
  	  'HWI-EAS255_8305_FC30G7J_PE:1:28:831:1267' => '1',
  	  'HWI-EAS255_8305_FC30G7J_PE:2:89:746:1917' => '1',
  	  'HWUSI-EAS493_8289_FC30GNW_PE:3:2:1019:1848' => '1',
  	  'HWI-EAS107_8284_FC30GCE_PE:1:66:685:1743' => '1',
  	  'USI-EAS39_8289_FC30GCV_PE:4:94:653:1757' => '1',
  	  'HWI-EAS255_8305_FC30G7J_PE:4:64:1773:283' => '1',
  	  'USI-EAS39_8289_FC30GCV_PE:6:57:1260:1513' => '1',
  	  'HWI-EAS255_8282_FC30G79_PE:6:49:937:1687' => '1',
  	  'HWI-EAS258_8310_FC30G9N_PE:5:21:240:1997' => '1',
  	  'USI-EAS39_8289_FC30GCV_PE:8:91:1422:1666' => '1',
  	  'HWI-EAS258_8282_FC30GAE_PE:2:8:439:239' => '1',
  	  'HWI-EAS107_8284_FC30GCE_PE:3:93:181:117' => '1',
  	  'HWI-EAS255_8282_FC30G79_PE:5:96:347:352' => '1',
  	  'HWI-EAS255_8282_FC30G79_PE:3:63:1555:951' => '1',
  	  'HWI-EAS300_8282_FC30BVC_PE:6:31:861:1658' => '1',
  	  'HWI-EAS255_8291_FC30GRN_PE:3:99:1600:752' => '1',
  	  'HWI-EAS258_8310_FC30G9N_PE:2:13:840:686' => '1',
  	  'HWI-EAS300_8282_FC30BVC_PE:1:67:664:610' => '1',
  	  'HWI-EAS255_8282_FC30G79_PE:3:95:1347:1361' => '1',
  	  'HWUSI-EAS493_8289_FC30GNW_PE:3:21:1321:1969' => '1',
  	  'HWI-EAS255_8305_FC30G7J_PE:6:61:1132:1881' => '1',
  	  'HWI-EAS255_8305_FC30G7J_PE:7:69:304:585' => '1',
  	  'USI-EAS39_8289_FC30GCV_PE:6:44:127:704' => '1',
  	},{
  	  'HWUSI-EAS493_8289_FC30GNW_PE:5:36:720:936' => '1',
  	  'HWI-EAS107_8284_FC30GCE_PE:6:48:1720:20' => '1',
  	  'HWI-EAS258_8310_FC30G9N_PE:2:58:1106:1805' => '1',
  	  'HWI-EAS138_4_FC30GP8:6:34:276:231' => '1',
  	  'HWI-EAS255_8282_FC30G79_PE:3:69:807:886' => '1',
  	  'HWI-EAS255_8282_FC30G79_PE:2:96:667:1546' => '1',
  	  'USI-EAS39_8289_FC30GCV_PE:6:46:858:1601' => '1',
  	  'HWI-EAS258_8310_FC30G9N_PE:1:100:66:1313' => '1',
  	  'HWUSI-EAS493_8289_FC30GNW_PE:4:88:71:146' => '1',
  	  'HWUSI-EAS493_8289_FC30GNW_PE:7:80:653:1227' => '1',
  	  'HWI-EAS258_8310_FC30G9N_PE:5:33:823:1300' => '1',
  	  'HWI-EAS255_8291_FC30GRN_PE:3:55:31:1737' => '1',
  	  'HWI-EAS255_8291_FC30GRN_PE:2:31:637:702' => '1',
  	  'HWI-EAS107_8284_FC30GCE_PE:5:87:1541:177' => '1',
  	  'HWUSI-EAS493_8289_FC30GNW_PE:2:82:186:1610' => '1',
  	  'HWI-EAS255_8282_FC30G79_PE:7:24:1637:969' => '1',
  	  'HWI-EAS138_4_FC30GP8:8:102:657:790' => '1',
  	  'HWI-EAS255_8291_FC30GRN_PE:5:95:923:1927' => '1',
  	  'HWI-EAS255_8291_FC30GRN_PE:8:85:705:1750' => '1',
  	  'HWI-EAS258_8310_FC30G9N_PE:8:92:1499:1751' => '1',
  	  'USI-EAS39_8289_FC30GCV_PE:3:38:1459:62' => '1',
  	  'HWI-EAS138_4_FC30GP8:3:33:332:188' => '1',
  	  'HWI-EAS255_8282_FC30G79_PE:3:37:624:1424' => '1',
  	  'HWI-EAS258_8310_FC30G9N_PE:2:42:1053:974' => '1',
  	  'HWI-EAS300_8282_FC30BVC_PE:6:78:1698:1141' => '1',
  	  'HWI-EAS255_8305_FC30G7J_PE:3:40:786:1068' => '1',
  	  'HWI-EAS255_8282_FC30G79_PE:2:10:1377:566' => '1',
  	  'HWI-EAS255_8305_FC30G7J_PE:3:19:1098:753' => '1',
  	  'HWI-EAS255_8291_FC30GRN_PE:5:87:1366:310' => '1',
  	  'HWI-EAS138_4_FC30GP8:4:19:1682:1439' => '1',
  	  'HWI-EAS255_8282_FC30G79_PE:7:42:818:255' => '1',
  	  'HWI-EAS300_8282_FC30BVC_PE:7:80:394:1321' => '1',
  	}];

  	is_deeply([Sanger::CGP::Pindel::OutputGen::CombinedRecordGenerator::_read_depth($wt_sam, 22, 16060480, 16060485)],$expect,'Checking correct interpretation of loci');
  };

  subtest '_count_sam_event_reads' => sub {

  	my $mt_sam = Bio::DB::HTS->new(-bam => $mt_file, -fasta => $ref_file);
  	my $samp_type_key = 'mt';

  	my $record_2 = new Sanger::CGP::Pindel::OutputGen::CombinedRecord(
  	  -length => 6,
      -type => 'D',
      -idx => 'D0',
      -alt_seq => undef,
      -chro => 22,
      -start => 16060480,
      -end => 16060485,
      -range_start => 16060479,
      -range_end => 16060525,
      -s1 => 12,
      -s2 => 234.764,
      -sum_ms => 202,
      -ref_seq => 'tttttc',
      -reads => {'test'=>{'-'=>[['EAS139_64:1:55:1728:1427_r1_D0',16,22,16060468,29,'12M6D63M','*',0,0,'AGTTAACTCTCTTTTTTCTTTTTCTTTTTCTTTTTCTTTTTCTTTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTT','*','MD:Z:1A10^TTTTTC63','NM:i:7'],]}},
      -lub => 'T'
  	);

    my $record_1 = new Sanger::CGP::Pindel::OutputGen::CombinedRecord(
      -length => 4,
      -type => 'I',
      -idx => 'I0',
      -alt_seq => 'AAAC',
      -chro => 22,
      -start => 16052167,
      -end => 16052167,
      -range_start => 16052167,
      -range_end => 16052200,
      -s1 => 16,
      -s2 => 550.185,
      -sum_ms => 533,
      -reads => {'test'=>{'-'=>[['EAS54_120:2:98:542:174_r1_I0',16,22,16052144,29,'12M6D63M','*',0,0,'ACAGAGCAAGACTCTATCTCAAAAAACAAACAAACAAACAAACAAACAAACAAACAAACTGTCAAAATCTGTAC','*','MD:Z:1A1063','NM:i:4'],]}},
      -lub => 'A'
    );

    my $exp_pos_reads = {
   	  'EAS25_5:1:71:437:248' => 1,
      'EAS192_65:5:10:134:1020' => 1,
      'EAS89_7:1:18:114:1055' => 1,
      'EAS89_7:2:24:1384:957' => 1,
      'EAS139_64:6:90:1172:620' => 1,
      'EAS131_8:3:53:1497:895' => 1,
      'EAS54_120:3:83:51:90' => 1,
      'EAS139_60:2:55:660:1032' => 1,
      'EAS54_120:1:14:1351:1949' => 1,
      'EAS192_60:3:72:1166:2006' => 1
    };

    my $exp_neg_reads = {
    	'EAS192_60:4:16:399:753' => 1,
        'EAS131_6:2:71:269:644' => 1,
        'EAS89_7:6:18:195:315' => 1
    };

  	my ($pos_reads,$neg_reads) = Sanger::CGP::Pindel::OutputGen::CombinedRecordGenerator::_count_sam_event_reads($record_1, $mt_sam, $samp_type_key);

  	is($record_1->b_mt_pos,10,'Check for correct pos read call counts');
  	is($record_1->b_mt_neg,3,'Check for correct neg read call counts');

  	is_deeply($pos_reads,$exp_pos_reads,'Check for correct pos read call hash');
  	is_deeply($neg_reads,$exp_neg_reads,'Check for correct neg read call hash');

  };

};

subtest 'Object funcions' => sub {

  subtest '_process_counts' => sub {


  	my $fai = Bio::DB::HTS::Fai->load($ref_file);
    my $wt_sam = Bio::DB::HTS->new(-bam => $wt_file, -fasta => $ref_file);
    my $mt_sam = Bio::DB::HTS->new(-bam => $mt_file, -fasta => $ref_file);
  	my $samp_type_key = 'mt';

  	my $generator = new Sanger::CGP::Pindel::OutputGen::CombinedRecordGenerator(
  	  -path => $pin_file,
  	  -fai => $fai,
  	  -wt_sam => $wt_sam,
  	  -mt_sam => $mt_sam,
  	  -mutant_sample_name => 'test'
  	);

  	my $record_1 = new Sanger::CGP::Pindel::OutputGen::CombinedRecord(
  	  -length => 6,
      -type => 'D',
      -idx => 'D0',
      -alt_seq => undef,
      -chro => 22,
      -start => 16060480,
      -end => 16060485,
      -range_start => 16060479,
      -range_end => 16060525,
      -ref_seq => 'tttttc',
      -reads => {
        'test'=>{'-'=>[['EAS139_64:1:55:1728:1427_r1_D0',16,22,16060468,29,'12M6D63M','*',0,0,'AGTTAACTCTCTTTTTTCTTTTTCTTTTTCTTTTTCTTTTTCTTTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTT','*','MD:Z:1A10^TTTTTC63','NM:i:7']]},
        'test-BL'=>{'-'=>[['EAS139_64:1:55:1728:1427_r1_D0',16,22,16060468,29,'12M6D63M','*',0,0,'AGTTAACTCTCTTTTTTCTTTTTCTTTTTCTTTTTCTTTTTCTTTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTT','*','MD:Z:1A10^TTTTTC63','NM:i:7']]}
      },
      -lub => 'T'
  	);

  	my $exp_record_1 = new Sanger::CGP::Pindel::OutputGen::CombinedRecord(
  	  -length => 6,
      -type => 'D',
      -idx => 'D0',
      -alt_seq => undef,
      -chro => 22,
      -start => 16060480,
      -end => 16060485,
      -range_start => 16060479,
      -range_end => 16060525,
      -ref_seq => 'tttttc',
      -reads => {
        'test'=>{'-'=>[['EAS139_64:1:55:1728:1427_r1_D0',16,22,16060468,29,'12M6D63M','*',0,0,'AGTTAACTCTCTTTTTTCTTTTTCTTTTTCTTTTTCTTTTTCTTTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTT','*','MD:Z:1A10^TTTTTC63','NM:i:7']]},
        'test-BL'=>{'-'=>[['EAS139_64:1:55:1728:1427_r1_D0',16,22,16060468,29,'12M6D63M','*',0,0,'AGTTAACTCTCTTTTTTCTTTTTCTTTTTCTTTTTCTTTTTCTTTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTT','*','MD:Z:1A10^TTTTTC63','NM:i:7']]}
      },
      -lub => 'T',

      -p_mt_pos => 0,
      -p_mt_neg => 1,
      -b_mt_pos => 0, ## no bwa calls which is a little odd..... as would expect there to be at least one bwa call in this case.....
      -b_mt_neg => 0, ## no bwa calls which is a little odd.....
      -d_mt_pos => 26,
      -d_mt_neg => 29,

      -rd_mt_pos => 26,
      -rd_mt_neg => 29,
      -uc_mt_pos => 0,
      -uc_mt_neg => 1,
      -fc_mt => 1,
      -fd_mt => 55,
  	);

  	$generator->_process_counts($record_1, $samp_type_key);

  	is_deeply($record_1,$exp_record_1,'Check for correct depth and call counts deletion');

  	my $record_2 = new Sanger::CGP::Pindel::OutputGen::CombinedRecord(
      -type => 'I',
      -idx => 'I1',
      -chro => 22,
      -start => 16052167,
      -end => 16052170,
      -range_start => 16052166,
      -range_end => 16052201,
      -ref_seq => undef,
      -alt_seq => 'AAAC',
      -reads => {
      	'test'=>{'+'=>[['EAS54_120:3:83:51:90_r1_I1',0,22,16052145,37,'23M4I47M','*',0,0,'ACAGAGCAAGACTCTATCTCAAAAAACAAACAAACAAACAAACAAACAAACAAACAAACTGTCAAAATCTGTAC','*','MD:Z:70','NM:i:4'],
      		           ['BOB_120:3:83:51:90_r1_I1',0,22,16052145,37,'23M4I47M','*',0,0,'ACAGAGCAAGACTCTATCTCAAAAAACAAACAAACAAACAAACAAACAAACAAACAAACTGTCAAAATCTGTAC','*','MD:Z:70','NM:i:4']]},
        'test-BL'=>{'+'=>[['EAS54_120:3:83:51:90_r1_I1',0,22,16052145,37,'23M4I47M','*',0,0,'ACAGAGCAAGACTCTATCTCAAAAAACAAACAAACAAACAAACAAACAAACAAACAAACTGTCAAAATCTGTAC','*','MD:Z:70','NM:i:4']]},
      }
    );

    my $exp_record_2 = new Sanger::CGP::Pindel::OutputGen::CombinedRecord(
      -type => 'I',
      -idx => 'I1',
      -chro => 22,
      -start => 16052167,
      -end => 16052170,
      -range_start => 16052166,
      -range_end => 16052201,
      -ref_seq => undef,
      -alt_seq => 'AAAC',
      -reads => {
      	'test'=>{'+'=>[['EAS54_120:3:83:51:90_r1_I1',0,22,16052145,37,'23M4I47M','*',0,0,'ACAGAGCAAGACTCTATCTCAAAAAACAAACAAACAAACAAACAAACAAACAAACAAACTGTCAAAATCTGTAC','*','MD:Z:70','NM:i:4'],
      		           ['BOB_120:3:83:51:90_r1_I1',0,22,16052145,37,'23M4I47M','*',0,0,'ACAGAGCAAGACTCTATCTCAAAAAACAAACAAACAAACAAACAAACAAACAAACAAACTGTCAAAATCTGTAC','*','MD:Z:70','NM:i:4']]},
        'test-BL'=>{'+'=>[['EAS54_120:3:83:51:90_r1_I1',0,22,16052145,37,'23M4I47M','*',0,0,'ACAGAGCAAGACTCTATCTCAAAAAACAAACAAACAAACAAACAAACAAACAAACAAACTGTCAAAATCTGTAC','*','MD:Z:70','NM:i:4']]}
      },

      -p_mt_pos => 2,
      -p_mt_neg => 0,
      -b_mt_pos => 10,
      -b_mt_neg => 3,
      -d_mt_pos => 34,
      -d_mt_neg => 28,

      -rd_mt_pos => 35,
      -rd_mt_neg => 28,
      -uc_mt_pos => 11, ## the pindel mapped read is presant in the mapped bam file but BOB is not...
      -uc_mt_neg => 3,
      -fc_mt => 14,
      -fd_mt => 63,
    );


    $generator->_process_counts($record_2, $samp_type_key);
    is_deeply($record_2,$exp_record_2,'Check for correct depth and call counts insertion');

  };


  subtest '_process_record' => sub {
    pass('WARNING! _process_record not actually tested, however individual components have been.');
  };
};

done_testing();
