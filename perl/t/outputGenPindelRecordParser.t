########## LICENCE ##########
# Copyright (c) 2014 Genome Research Ltd. 
#  
# Author: Keiran Raine <cgpit@sanger.ac.uk> 
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

use Bio::DB::Sam;

my $MODULE = 'Sanger::CGP::Pindel::OutputGen::PindelRecordParser';

use_ok 'Sanger::CGP::Vcf::Sample';
use_ok 'Sanger::CGP::Vcf::Contig';
use_ok 'Sanger::CGP::Pindel::OutputGen::PindelRecord';

my $test_data = "$Bin/data";
my $ref_file = File::Spec->catfile($test_data,'genome_22.fa');
my $pin_file = File::Spec->catfile($test_data,'test.txt_22');

subtest 'Initialisation checks' => sub {
  use_ok($MODULE);

  my $fai = Bio::DB::Sam::Fai->load($ref_file);
  my $obj = new_ok($MODULE => [-path => $pin_file, -fai => $fai]); ## this will fail without a valid sam pindel file

  open(my $fh,"<",$pin_file);
  $obj = new_ok($MODULE => [-fh => $fh, -fai => $fai]); ## this will fail without a valid sam pindel file
  close($fh);

  is($obj->{_fh},$fh,'Correctly setting fh object');
  is($obj->{_fai},$fai,'Correctly setting fai object');

  eval{
    new $MODULE(-fai => $fai);
  };unless($@){
  	fail('Initialisation checks: expected an exception to be thrown if no pindel file path set.');
  }else{
  	pass('new method throws exception as expected if no -path or -fh arguments');
  }
};

subtest 'Non-object funcions' => sub {

  subtest '_shrink_change' => sub {

    my $max_repeat_unit = 100;
    my $test1 = 'ATATATATATATATATA';
    my $exp1 = 'ATATATATATATATATA';
    is(Sanger::CGP::Pindel::OutputGen::PindelRecordParser::_shrink_change($test1,$max_repeat_unit),$exp1,'Unbalanced AT repeat');

    my $test2 = 'ATATATGTATATATAT';
    my $exp2 = 'ATATATGTATATATAT';
    is(Sanger::CGP::Pindel::OutputGen::PindelRecordParser::_shrink_change($test2,$max_repeat_unit),$exp2,'Imperfect AT repeat');

    my $test3 = 'AAAAAAA';
    my $exp3 = 'A';
    is(Sanger::CGP::Pindel::OutputGen::PindelRecordParser::_shrink_change($test3,$max_repeat_unit),$exp3,'Polly A');

    my $test4 = 'ATATATATATATATAT';
    my $exp4 = 'AT';
    is(Sanger::CGP::Pindel::OutputGen::PindelRecordParser::_shrink_change($test4,$max_repeat_unit),$exp4,'AT repeat');

    my $test5 = 'ATCATCATCATC';
    my $exp5 = 'ATC';
    is(Sanger::CGP::Pindel::OutputGen::PindelRecordParser::_shrink_change($test5,$max_repeat_unit),$exp5,'ATC repeat');

    my $test6 = 'ATCATGATCATG';
    my $exp6 = 'ATCATG';
    is(Sanger::CGP::Pindel::OutputGen::PindelRecordParser::_shrink_change($test6,$max_repeat_unit),$exp6,'ATCATG dupe');

    $max_repeat_unit = 3;

    my $test7 = 'ATCATCATC';
    my $exp7 = 'ATC';
    is(Sanger::CGP::Pindel::OutputGen::PindelRecordParser::_shrink_change($test7,$max_repeat_unit),$exp7,'ATC repeat with max_repeat_unit at 3');

    my $test8 = 'ATCTATCT';
    my $exp8 = 'ATCTATCT';
    is(Sanger::CGP::Pindel::OutputGen::PindelRecordParser::_shrink_change($test8,$max_repeat_unit),$exp8,'ATCT repeat with max_repeat_unit at 3');
  };

  subtest '_parse_header' => sub {

    my $header_string_1 = q(0	D 6	NT 0 ""	ChrID 22	BP 16060479	16060486	BP_range 16060479	16060525	Supports 6	+ 5	- 1	S1 12	S2 234.764	SUM_MS 202	NumSupSamples 2	test 4	test-BL 2);
    my $record_1 = new Sanger::CGP::Pindel::OutputGen::PindelRecord();
    my $exp_record_1 = new Sanger::CGP::Pindel::OutputGen::PindelRecord(
      -length => 6,
      -type => 'D',
      -idx => 'D0',
      -alt_seq => undef,
      -chro => 22,
      -start => 16060480,
      -end => 16060485,
      -range_start => 16060479,
      -range_end => 16060525,
      -num_samples => 2,
      -sample_contrib => {'test'    => '4',
                          'test-BL' => '2'},
      -s1 => 12,
      -s2 => 234.764,
      -sum_ms => 202
    );
    is_deeply(Sanger::CGP::Pindel::OutputGen::PindelRecordParser::_parse_header($record_1,$header_string_1),$exp_record_1,'Deletion variant header');

    my $header_string_2 = q(0	I 4	NT 4 "AAAC"	ChrID 22	BP 16052167	16052168	BP_range 16052167	16052200	Supports 15	+ 0	- 15	S1 16	S2 550.185	SUM_MS 533	NumSupSamples 2	test 9	test-BL 6);
    my $record_2 = new Sanger::CGP::Pindel::OutputGen::PindelRecord();
    my $exp_record_2 = new Sanger::CGP::Pindel::OutputGen::PindelRecord(
      -length => 4,
      -type => 'I',
      -idx => 'I0',
      -alt_seq => 'AAAC',
      -chro => 22,
      -start => 16052167,
      -end => 16052167,
      -range_start => 16052167,
      -range_end => 16052200,
      -num_samples => 2,
      -sample_contrib => {'test'    => '9',
                          'test-BL' => '6'},
      -s1 => 16,
      -s2 => 550.185,
      -sum_ms => 533
    );
    is_deeply(Sanger::CGP::Pindel::OutputGen::PindelRecordParser::_parse_header($record_2,$header_string_2),$exp_record_2,'Insertion variant header');

    my $header_string_3 = q(8161	D 18	NT 12 "GCCCCCCATCGC"	ChrID 22	BP 51154501	51154520	Supports 3	+ 0	- 3	S1 4	SUM_MS 180	NumSupSamples 2	test 2	test-BL 1);
    my $record_3 = new Sanger::CGP::Pindel::OutputGen::PindelRecord();
    my $exp_record_3 = new Sanger::CGP::Pindel::OutputGen::PindelRecord(
      -length => 18,
      -type => 'DI',
      -idx => 'DI8161',
      -alt_seq => 'GCCCCCCATCGC',
      -chro => 22,
      -start => 51154502,
      -end => 51154519,
      -range_start => 51154501,
      -range_end => 51154520,
      -num_samples => 2,
      -sample_contrib => {'test'    => '2',
                          'test-BL' => '1'},
      -s1 => 4,
      -sum_ms => 180
    );
    is_deeply(Sanger::CGP::Pindel::OutputGen::PindelRecordParser::_parse_header($record_3,$header_string_3),$exp_record_3,'Complex variant header');

  };

  subtest '_parse_header_v02' => sub {

    my $header_string_1 = q(244	D 1671	NT 0 ""	ChrID 20	BP 1337143	1338815	BP_range 1337143	1338818	Supports 10	10	+ 6	6	- 4	4	S1 35	SUM_MS 322	2	NumSupSamples 2 2	COLO-829 3 3 0 0	COLO-829-BL 3 3 4 4
);
    my $record_1 = new Sanger::CGP::Pindel::OutputGen::PindelRecord();
    my $exp_record_1 = new Sanger::CGP::Pindel::OutputGen::PindelRecord(
      -length => 1671,
      -type => 'D',
      -idx => 'D244',
      -alt_seq => undef,
      -chro => 20,
      -start => 1337144,
      -end => 1338814,
      -range_start => 1337143,
      -range_end => 1338818,
      -s1 => 35,
      -sum_ms => 322
    );
    is_deeply(Sanger::CGP::Pindel::OutputGen::PindelRecordParser::_parse_header_v02($record_1,$header_string_1),$exp_record_1,'Deletion variant header');

    my $header_string_2 = q(57	I 10	NT 10 "ACTTGTTCCC"	ChrID 20	BP 400747	400748	BP_range 400742	400749	Supports 3	2	+ 0	0	- 3	2	S1 4	SUM_MS 111	2	NumSupSamples 2	2	COLO-829 0 0 2 1	COLO-829-BL 0 0 1 1);
    my $record_2 = new Sanger::CGP::Pindel::OutputGen::PindelRecord();
    my $exp_record_2 = new Sanger::CGP::Pindel::OutputGen::PindelRecord(
      -length => 10,
      -type => 'I',
      -idx => 'I57',
      -alt_seq => 'ACTTGTTCCC',
      -chro => 20,
      -start => 400747,
      -end => 400747,
      -range_start => 400742,
      -range_end => 400749,
      -s1 => 4,
      -sum_ms => 111
    );
    is_deeply(Sanger::CGP::Pindel::OutputGen::PindelRecordParser::_parse_header_v02($record_2,$header_string_2),$exp_record_2,'Insertion variant header');

    my $header_string_3 = q(57	D 10	NT 10 "ACTTGTTCCC"	ChrID 20	BP 400747	400758	Supports 3	2	+ 0	0	- 3	2	S1 4	SUM_MS 111	2	NumSupSamples 2	2	COLO-829 0 0 2 1	COLO-829-BL 0 0 1 1);
    my $record_3 = new Sanger::CGP::Pindel::OutputGen::PindelRecord();
    my $exp_record_3 = new Sanger::CGP::Pindel::OutputGen::PindelRecord(
      -length => 10,
      -type => 'DI',
      -idx => 'DI57',
      -alt_seq => 'ACTTGTTCCC',
      -chro => 20,
      -start => 400748,
      -end => 400757,
      -range_start => 400748,
      -range_end => 400757,
      -s1 => 4,
      -sum_ms => 111
    );
    is_deeply(Sanger::CGP::Pindel::OutputGen::PindelRecordParser::_parse_header_v02($record_3,$header_string_3),$exp_record_3,'Complex variant header');

  };

  subtest '_repeat_count' => sub {

    my $left = 'ATCGTGTGTCCGCAAAAA';
    my $right = 'TCTCAAAAAAAAAAAAAAAAAA';

    my $record_1 = new Sanger::CGP::Pindel::OutputGen::PindelRecord(
      -start => 16060480,
      -end => 16060484,
      -range_start => 16060479,
      -range_end => 16060489,
      -ref_seq =>'TCTC',
      -alt_seq => undef,
      -min_change => 'TC',
  	);

    is_deeply(Sanger::CGP::Pindel::OutputGen::PindelRecordParser::_repeat_count($record_1,\$left,\$right),4,'Check deletion repeats');

    my $record_2 = new Sanger::CGP::Pindel::OutputGen::PindelRecord(
      -start => 16060480,
      -end => 16060484,
      -range_start => 16060479,
      -range_end => 16060489,
      -ref_seq => undef,
      -alt_seq => 'TCTC',
      -min_change => 'TC',
  	);

    is_deeply(Sanger::CGP::Pindel::OutputGen::PindelRecordParser::_repeat_count($record_2,\$left,\$right),2,'Check insertion repeats');

    my $record_3 = new Sanger::CGP::Pindel::OutputGen::PindelRecord(
      -start => 16060480,
      -end => 16060484,
      -range_start => 16060479,
      -range_end => 16060489,
      -ref_seq => 'TGTG',
      -alt_seq => 'TCTC',
      -min_change => 'TC',
  	);

    is_deeply(Sanger::CGP::Pindel::OutputGen::PindelRecordParser::_repeat_count($record_3,\$left,\$right),2,'Check complex repeats');

  };

  subtest 'calmd' => sub {

    my $fai = Bio::DB::Sam::Fai->load($ref_file);

    ##ÊThese are real coordinates from human build 37
    #genome.fa 22:30000000-30000050
    #ATCGCTTCCCGCATGAGCTTCAGCTCTCTCAAGAGGAAGCAACCCAAGACG

    my $chr = 22;
    my $start = 30000002;
    #my $cigar1 = '2S20M3I25M6S';
    my $cigar1 = [qw{2 S 20 M 3 I 25 M 6 S}];
    my $seq_ref1 = 'ATCGCTTCCCGCATGAGCTTCATTTGCTCTCTCAAGAGGAAGCAACCCAAGGGGGG';
    my $ref_sequence = $fai->fetch($chr.':30000000-30001000');
    my $ref_start = 30000000;

    my $exp_passed_to_fai_arg_value1 = []; ## insertion nothing should be passed...
    my @exp1 = qw{MD:Z:45 NM:i:3 30000046};

    is_deeply([Sanger::CGP::Pindel::OutputGen::PindelRecordParser::calmd($chr, $start, $cigar1, \$seq_ref1,\$ref_sequence,$ref_start)],\@exp1,'calmd insertion');

    #my $cigar2 = '2S20M3D25M6S';
    my $cigar2 = [qw{2 S 20 M 3 D 25 M 6 S}];
    my $seq_ref2 = 'ATCGCTTCCCGCATGAGCTTCACTCTCAAGAGGAAGCAACCCAAGACGGGGGG';
    my @exp2 = qw{MD:Z:20^GCT25 NM:i:3 30000049};

    is_deeply([Sanger::CGP::Pindel::OutputGen::PindelRecordParser::calmd($chr, $start, $cigar2, \$seq_ref2,\$ref_sequence,$ref_start)],\@exp2,'calmd deletion');

    #my $cigar3 = '2S20M3D25M6S';
    my $cigar3 = [qw{2 S 20 M 3 D 25 M 6 S}];
    my $seq_ref3 = 'ATCGaTTCCCGCATGAGCTTCACTCTCAAGAGGAAGCAACCCAAGACGGGGGG';
    my @exp3 = qw{MD:Z:2C17^GCT25 NM:i:4 30000049};

    is_deeply([Sanger::CGP::Pindel::OutputGen::PindelRecordParser::calmd($chr, $start, $cigar3, \$seq_ref3,\$ref_sequence,$ref_start)],\@exp3,'calmd sub deletion');

  };

  subtest '_parse_read' => sub {

    my $fai = Bio::DB::Sam::Fai->load($ref_file);

    my $ref_sequence_a = $fai->fetch('22:16060000-16061000');
    my $ref_start_a = 16060000;

    my $ref_sequence_b = $fai->fetch('22:16052000-16053000');
    my $ref_start_b = 16052000;

    my $ref_sequence_c = $fai->fetch('22:51154000-51155000');
    my $ref_start_c = 51154000;

    ## Actually contains a sub and a deletion of tttttc
    my $read_1 = q(                                                               AGTTAACTCTCT      TTTTTCTTTTTCTTTTTCTTTTTCTTTTTCTTTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTT            	+	16060230	29	test	@EAS139_64:1:55:1728:1427/2);
    my $change_ref_start_1 = 75;
    my $change_ref_end_1 = 81;
    my $record_1 = new Sanger::CGP::Pindel::OutputGen::PindelRecord(
      -type => 'D',
      -idx => 'D1',
      -chro => 22,
      -start => 16060480,
      -ref_seq => 'tttttc'
    );

    my $exp_record1 = new Sanger::CGP::Pindel::OutputGen::PindelRecord(
      -type => 'D',
      -idx => 'D1',
      -chro => 22,
      -start => 16060480,
      -ref_seq => 'tttttc',
      -reads => {'test'=>{'-'=>[['EAS139_64:1:55:1728:1427_r1_D1',16,22,16060468,29,'12M6D63M','*',0,0,'AGTTAACTCTCTTTTTTCTTTTTCTTTTTCTTTTTCTTTTTCTTTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTT','*','MD:Z:1A10^TTTTTC63','NM:i:7']]}}
    );
    # $name,$flag,$record->chro(),$start_pos,$mapq,$cigar,'*','0','0',$read_seq,'*',$read_group,@tags

    Sanger::CGP::Pindel::OutputGen::PindelRecordParser::_parse_read($record_1, $record_1->chro, $record_1->start, \$read_1, length($record_1->ref_seq),('1_D1'), $change_ref_start_1, $change_ref_end_1,\$ref_sequence_a,$ref_start_a);
    is_deeply($record_1,$exp_record1,'Read with a small deletion');

   ## Actually contains a sub and a deletion of tttttc
    my $read_2 = q(                              TCTTTTTTCTTTTTCTTTTTCTTTTTCTTTTTCTTTTTCTTTTTC              TCTTTCTTTCTTTCTTTCTTTCTTTCTTTC                                             	+	16060166	60	test	@EAS131_6:8:80:742:1825/1);
    my $change_ref_start_2 = 75;
    my $change_ref_end_2 = 89;
    my $record_2 = new Sanger::CGP::Pindel::OutputGen::PindelRecord(
      -type => 'D',
      -idx => 'D1',
      -chro => 22,
      -start => 16060522,
      -ref_seq => 'TTTCTTTCTTTCTTTCTTTCTTTCTTTCTT'
    );

    ## Actually contains a large deletion and a sub right at the end...
    my $exp_record2 = new Sanger::CGP::Pindel::OutputGen::PindelRecord(
      -type => 'D',
      -idx => 'D1',
      -chro => 22,
      -start => 16060522,
      -ref_seq => 'TTTCTTTCTTTCTTTCTTTCTTTCTTTCTT',
      -reads => {'test'=>{'-'=>[['EAS131_6:8:80:742:1825_r1_D1',16,22,16060477,60,'45M30D30M','*',0,0,'TCTTTTTTCTTTTTCTTTTTCTTTTTCTTTTTCTTTTTCTTTTTCTCTTTCTTTCTTTCTTTCTTTCTTTCTTTC','*','MD:Z:45^TTTCTTTCTTTCTTTCTTTCTTTCTTTCTT29T0','NM:i:31']]}}
    );
    # $name,$flag,$record->chro(),$start_pos,$mapq,$cigar,'*','0','0',$read_seq,'*',$read_group,@tags

    Sanger::CGP::Pindel::OutputGen::PindelRecordParser::_parse_read($record_2, $record_2->chro, $record_2->start, \$read_2, length($record_2->ref_seq),('1_D1'), $change_ref_start_2, $change_ref_end_2,\$ref_sequence_a,$ref_start_a);
    is_deeply($record_2,$exp_record2,'Read with a large deletion');

    ## Actually contains an insertion
    my $read_3 = q(                                                    ACAGAGCAAGACTCTATCTCAAAAAACAAACAAACAAACAAACAAACAAACAAACAAACTGTCAAAATCTGTAC                        	-	16052456	37	test	@EAS54_120:3:83:51:90/2
);
    my $change_ref_start_3 = 75;
    my $change_ref_end_3 = 79;
    my $record_3 = new Sanger::CGP::Pindel::OutputGen::PindelRecord(
      -type => 'I',
      -idx => 'I1',
      -chro => 22,
      -start => 16052167,
      -alt_seq => 'AAAC'
    );

    my $exp_record3 = new Sanger::CGP::Pindel::OutputGen::PindelRecord(
      -type => 'I',
      -idx => 'I1',
      -chro => 22,
      -start => 16052167,
      -alt_seq => 'AAAC',
      -reads => {'test'=>{'+'=>[['EAS54_120:3:83:51:90_r1_I1',0,22,16052145,37,'23M4I47M','*',0,0,'ACAGAGCAAGACTCTATCTCAAAAAACAAACAAACAAACAAACAAACAAACAAACAAACTGTCAAAATCTGTAC','*','MD:Z:70','NM:i:4']]}}
    );
    # $name,$flag,$record->chro(),$start_pos,$mapq,$cigar,'*','0','0',$read_seq,'*',$read_group,@tags

    Sanger::CGP::Pindel::OutputGen::PindelRecordParser::_parse_read($record_3, $record_3->chro, ($record_3->start + 1), \$read_3, length($record_3->ref_seq),('1_I1'), $change_ref_start_3, $change_ref_end_3,\$ref_sequence_b,$ref_start_b);
    is_deeply($record_3,$exp_record3,'Read with an insertion');

    ## Actually contains a complex variant with substitutuions at either end
    my $read_4 = q(                                                           CTGTCCAGTGGGCACCGCCCCCCATCGCCTCATCCCTCCCATGGGCAGTCTCATCCCTCTCCCCAGCTGCCCCTC		-	51154758	60	test	@EAS139_60:6:82:717:1216/2
);
    my $change_ref_start_4 = 75;
    my $change_ref_end_4 = 87;
    my $record_4 = new Sanger::CGP::Pindel::OutputGen::PindelRecord(
      -type => 'DI',
      -idx => 'DI1',
      -chro => 22,
      -start => 51154502,
      -ref_seq => 'CCCACCTCCCCATCACCT',
      -alt_seq => 'GCCCCCCATCGC'
    );

    my $exp_record4 = new Sanger::CGP::Pindel::OutputGen::PindelRecord(
      -type => 'DI',
      -idx => 'DI1',
      -chro => 22,
      -start => 51154502,
      -ref_seq => 'CCCACCTCCCCATCACCT',
      -alt_seq => 'GCCCCCCATCGC',
      -reads => {'test'=>{'+'=>[['EAS139_60:6:82:717:1216_r1_DI1',0,22,51154486,60,'16M18D12I47M','*',0,0,'CTGTCCAGTGGGCACCGCCCCCCATCGCCTCATCCCTCCCATGGGCAGTCTCATCCCTCTCCCCAGCTGCCCCTC','*','MD:Z:7C8^CCCACCTCCCCATCACCT20C9G12A3','NM:i:34']]}}
    );
    # $name,$flag,$record->chro(),$start_pos,$mapq,$cigar,'*','0','0',$read_seq,'*',$read_group,@tags

    Sanger::CGP::Pindel::OutputGen::PindelRecordParser::_parse_read($record_4, $record_4->chro, $record_4->start, \$read_4, length($record_4->ref_seq),('1_DI1'), $change_ref_start_4, $change_ref_end_4,\$ref_sequence_c,$ref_start_c);
    is_deeply($record_4,$exp_record4,'Read with a complex');
  };

};

subtest 'Object funcions' => sub {

  subtest '_parse_aignment' => sub {

  	my $fai = Bio::DB::Sam::Fai->load($ref_file);
    my $obj = new_ok($MODULE => [-path => $pin_file, -fai => $fai]); ## this will fail without a valid sam pindel file

  	my $header_string_1 = q(GAGACCTCCCCAGAAATGGATGCCAGCATTATGCTTCCTATACAGCCTGCAGAACCATGAGCCAATTAACTCTCTtttttcTTTTTCTTTTTCTTTTTCTTTTTCTTTTTCTTTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCT);
    my $alignments_1 = [
         q{                                                               AGTTAACTCTCT      TTTTTCTTTTTCTTTTTCTTTTTCTTTTTCTTTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTT            	+	16060230	29	test	@EAS139_64:1:55:1728:1427/2},
         q{                                                               AGTTAACTCTCT      TTTTTCTTTTTCTTTTTCTTTTTCTTTTTCTTTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTT            	+	16060230	29	test	@EAS139_64:1:55:1728:1427/2},
       ];

  	my $record_1 = new Sanger::CGP::Pindel::OutputGen::PindelRecord(
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
      -ref_seq => lc'TTTTTC',
  	);

    my $exp_record_1 = new Sanger::CGP::Pindel::OutputGen::PindelRecord(
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
      -ref_seq => lc'TTTTTC',
      -reads => {'test'=>{'-'=>[['EAS139_64:1:55:1728:1427_r1_D0',16,22,16060468,29,'12M6D63M','*',0,0,'AGTTAACTCTCTTTTTTCTTTTTCTTTTTCTTTTTCTTTTTCTTTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTT','*','MD:Z:1A10^TTTTTC63','NM:i:7'],
      	                        ['EAS139_64:1:55:1728:1427_r2_D0',16,22,16060468,29,'12M6D63M','*',0,0,'AGTTAACTCTCTTTTTTCTTTTTCTTTTTCTTTTTCTTTTTCTTTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTT','*','MD:Z:1A10^TTTTTC63','NM:i:7']]}},
      -lub => 'T',
      -min_change => lc'TTTTTC',
      -repeats => 7
    );

    $obj->_parse_alignment($record_1,$alignments_1,\$header_string_1);
    is_deeply($record_1,$exp_record_1,'Check full small deletion record');

    my $header_string_2 = q(CAGCCTGCAGAACCATGAGCCAATTAACTCTCTTTTTTCTTTTTCTTTTTCTTTTTCTTTTTCTTTTTCTTTTTCtttct<20>ttcttTCTTTCTTTCTTTCTTTCTTTCTTTCTTTTTGTTTTCTTTCATCTTTCCTTCTTCTTTTTT);
    my $alignments_2 = [
         q{                              TCTTTTTTCTTTTTCTTTTTCTTTTTCTTTTTCTTTTTCTTTTTC              TCTTTCTTTCTTTCTTTCTTTCTTTCTTTC                                             	+	16060166	60	test	@EAS131_6:8:80:742:1825/1},
         q{                              TCTTTTTTCTTTTTCTTTTTCTTTTTCTTTTTCTTTTTCTTTTTC              TCTTTCTTTCTTTCTTTCTTTCTTTCTTTC                                             	+	16060166	60	test	@EAS131_6:8:80:742:1825/1},
       ];

  	my $record_2 = new Sanger::CGP::Pindel::OutputGen::PindelRecord(
  	  -length => 6,
  	  -type => 'D',
      -idx => 'D0',
      -chro => 22,
      -start => 16060522,
      -end => 16060551,
      -range_start => 16060521,
      -range_end => 16060553,
      -s1 => 12,
      -s2 => 234.764,
      -sum_ms => 202,
  	);

    my $exp_record_2 = new Sanger::CGP::Pindel::OutputGen::PindelRecord(
      -length => 6,
      -type => 'D',
      -idx => 'D0',
      -alt_seq => undef,
      -chro => 22,
      -start => 16060522,
      -end => 16060551,
      -range_start => 16060521,
      -range_end => 16060553,
      -s1 => 12,
      -s2 => 234.764,
      -sum_ms => 202,
      -ref_seq => 'TTTCTTTCTTTCTTTCTTTCTTTCTTTCTT',
      -reads => {'test'=>{'-'=>[['EAS131_6:8:80:742:1825_r1_D0',16,22,16060477,60,'45M30D30M','*',0,0,'TCTTTTTTCTTTTTCTTTTTCTTTTTCTTTTTCTTTTTCTTTTTCTCTTTCTTTCTTTCTTTCTTTCTTTCTTTC','*','MD:Z:45^TTTCTTTCTTTCTTTCTTTCTTTCTTTCTT29T0','NM:i:31'],
      	                        ['EAS131_6:8:80:742:1825_r2_D0',16,22,16060477,60,'45M30D30M','*',0,0,'TCTTTTTTCTTTTTCTTTTTCTTTTTCTTTTTCTTTTTCTTTTTCTCTTTCTTTCTTTCTTTCTTTCTTTCTTTC','*','MD:Z:45^TTTCTTTCTTTCTTTCTTTCTTTCTTTCTT29T0','NM:i:31']]}},
      -lub => 'C',
      -min_change => 'TTTCTTTCTTTCTTTCTTTCTTTCTTTCTT',
      -repeats => 1 # The minimum repeat is not actually repeated within the local vicinity of the event, despite the region being a repeat-region.
    );

    $obj->_parse_alignment($record_2,$alignments_2,\$header_string_2);
    is_deeply($record_2,$exp_record_2,'Check full large deletion record');

  };

  subtest '_process_record' => sub {
    pass('WARNING! _process_record not actually tested, however individual components have been.');
  };
};

done_testing();
