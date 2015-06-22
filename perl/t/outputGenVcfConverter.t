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
use Data::Dumper;

my $MODULE = 'Sanger::CGP::Pindel::OutputGen::VcfConverter';

use DateTime qw();

use_ok 'Sanger::CGP::Vcf::Sample';
use_ok 'Sanger::CGP::Vcf::Contig';
use_ok 'Sanger::CGP::Vcf::VcfProcessLog';
use_ok 'Sanger::CGP::Vcf::VcfUtil';
use_ok 'Sanger::CGP::Pindel::OutputGen::CombinedRecord';

subtest 'Initialisation checks' => sub {
  use_ok($MODULE);
  my $obj = new_ok($MODULE => []); ## this will fail without a valid sam file/object

  $obj = new_ok($MODULE => []); ## this will fail without a valid sam file/object
};

subtest 'Object funcions' => sub {

  subtest 'gen_header' => sub {

  	my $wt_sample = new Sanger::CGP::Vcf::Sample(-name => 'test_wt');
  	my $mt_sample = new Sanger::CGP::Vcf::Sample(-name => 'test_mt');

  	my $contig_1 = new Sanger::CGP::Vcf::Contig(-name => '1', -length => 10, -species => 'jelly fish', -assembly => '25');
  	my $contig_x = new Sanger::CGP::Vcf::Contig(-name => 'x', -length => 10, -species => 'jelly fish', -assembly => '25');

  	my @process_logs = ();

    push @process_logs, new Sanger::CGP::Vcf::VcfProcessLog(
	  -input_vcf_source => 'Pindel',
      -input_vcf_ver => 'test_version1',
	);

	push @process_logs, new Sanger::CGP::Vcf::VcfProcessLog(
	  -input_vcf_source => 'Jeff K',
      -input_vcf_ver => 'test_version2',
      -input_vcf_params => {'o'=>'my_file'},
	);

  	my $contigs = [$contig_1,$contig_x];

  	my $converter = new Sanger::CGP::Pindel::OutputGen::VcfConverter(-contigs => $contigs);

  	my $reference_name = 'BOB_ref';
    my $input_source = 'BOB';
  	my $date = DateTime->now->strftime('%Y%m%d');

  	my $exp1 =
qq{##fileformat=VCFv4.1
##fileDate=$date
##source_$date.1=BOB
##reference=BOB_ref
##contig=<ID=1,assembly=25,length=10,species="jelly fish">
##contig=<ID=x,assembly=25,length=10,species="jelly fish">
##INFO=<ID=PC,Number=1,Type=String,Description="Pindel call">
##INFO=<ID=RS,Number=1,Type=Integer,Description="Range start">
##INFO=<ID=RE,Number=1,Type=Integer,Description="Range end">
##INFO=<ID=LEN,Number=1,Type=Integer,Description="Length">
##INFO=<ID=REP,Number=1,Type=Integer,Description="Change repeat count within range">
##INFO=<ID=S1,Number=1,Type=Integer,Description="S1">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PP,Number=1,Type=Integer,Description="Pindel calls on the positive strand">
##FORMAT=<ID=NP,Number=1,Type=Integer,Description="Pindel calls on the negative strand">
##FORMAT=<ID=PB,Number=1,Type=Integer,Description="BWA calls on the positive strand">
##FORMAT=<ID=NB,Number=1,Type=Integer,Description="BWA calls on the negative strand">
##FORMAT=<ID=PD,Number=1,Type=Integer,Description="BWA mapped reads on the positive strand">
##FORMAT=<ID=ND,Number=1,Type=Integer,Description="BWA mapped reads on the negative strand">
##FORMAT=<ID=PR,Number=1,Type=Integer,Description="Total mapped reads on the positive strand">
##FORMAT=<ID=NR,Number=1,Type=Integer,Description="Total mapped reads on the negative strand">
##FORMAT=<ID=PU,Number=1,Type=Integer,Description="Unique calls on the positive strand">
##FORMAT=<ID=NU,Number=1,Type=Integer,Description="Unique calls on the negative strand">
##FORMAT=<ID=TG,Number=1,Type=Integer,Description="Total distinct contributing read groups">
##FORMAT=<ID=VG,Number=1,Type=Integer,Description="Variant distinct contributing read groups">
##vcfProcessLog=<InputVCFSource=<Pindel>,InputVCFVer=<test_version1>>
##vcfProcessLog=<InputVCFSource=<Jeff K>,InputVCFVer=<test_version2>,InputVCFParam=<o=my_file>>
##SAMPLE=<ID=NORMAL,SampleName=test_wt>
##SAMPLE=<ID=TUMOUR,SampleName=test_mt>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NORMAL	TUMOUR
};

  	my $exp2 =
qq{##fileformat=VCFv4.1
##fileDate=$date
##source_$date.1=BOB
##reference=BOB_ref
##contig=<ID=1,assembly=25,length=10,species="jelly fish">
##contig=<ID=x,assembly=25,length=10,species="jelly fish">
##INFO=<ID=PC,Number=1,Type=String,Description="Pindel call">
##INFO=<ID=RS,Number=1,Type=Integer,Description="Range start">
##INFO=<ID=RE,Number=1,Type=Integer,Description="Range end">
##INFO=<ID=LEN,Number=1,Type=Integer,Description="Length">
##INFO=<ID=REP,Number=1,Type=Integer,Description="Change repeat count within range">
##INFO=<ID=S1,Number=1,Type=Integer,Description="S1">
##INFO=<ID=S2,Number=1,Type=Float,Description="S2">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PP,Number=1,Type=Integer,Description="Pindel calls on the positive strand">
##FORMAT=<ID=NP,Number=1,Type=Integer,Description="Pindel calls on the negative strand">
##FORMAT=<ID=PB,Number=1,Type=Integer,Description="BWA calls on the positive strand">
##FORMAT=<ID=NB,Number=1,Type=Integer,Description="BWA calls on the negative strand">
##FORMAT=<ID=PD,Number=1,Type=Integer,Description="BWA mapped reads on the positive strand">
##FORMAT=<ID=ND,Number=1,Type=Integer,Description="BWA mapped reads on the negative strand">
##FORMAT=<ID=PR,Number=1,Type=Integer,Description="Total mapped reads on the positive strand">
##FORMAT=<ID=NR,Number=1,Type=Integer,Description="Total mapped reads on the negative strand">
##FORMAT=<ID=PU,Number=1,Type=Integer,Description="Unique calls on the positive strand">
##FORMAT=<ID=NU,Number=1,Type=Integer,Description="Unique calls on the negative strand">
##FORMAT=<ID=TG,Number=1,Type=Integer,Description="Total distinct contributing read groups">
##FORMAT=<ID=VG,Number=1,Type=Integer,Description="Variant distinct contributing read groups">
##vcfProcessLog=<InputVCFSource=<Pindel>,InputVCFVer=<test_version1>>
##vcfProcessLog=<InputVCFSource=<Jeff K>,InputVCFVer=<test_version2>,InputVCFParam=<o=my_file>>
##SAMPLE=<ID=NORMAL,SampleName=test_wt>
##SAMPLE=<ID=TUMOUR,SampleName=test_mt>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NORMAL	TUMOUR
};


    my $header_string_no_s2 = $converter->gen_header($wt_sample, $mt_sample, \@process_logs, 0, $reference_name, $input_source);
    my @header_string_no_s2_bits = split "\n", $header_string_no_s2;
   	my $header_string_no_s2_columns = pop @header_string_no_s2_bits;
   	my @exp1_bits = split "\n", $exp1;
   	my $exp1_columns = pop @exp1_bits;

   	is($header_string_no_s2_columns,$exp1_columns,'Check no S2 columns line');
   	is(scalar @header_string_no_s2_bits, scalar @exp1_bits,'Check no S2 correct number of lines');

   	foreach my $got_line (@header_string_no_s2_bits){
   		my $expected_line = shift @exp1_bits;
   		is_deeply(header_to_hash($got_line),header_to_hash($expected_line),'Checking no S2 header line');
   	}

    my $header_string_s2 = $converter->gen_header($wt_sample, $mt_sample, \@process_logs, 1, $reference_name, $input_source);

    my @header_string_s2_bits = split "\n", $header_string_s2;
   	my $header_string_s2_columns = pop @header_string_s2_bits;
   	my @exp2_bits = split "\n", $exp2;
   	my $exp2_columns = pop @exp2_bits;

   	is($header_string_s2_columns,$exp2_columns,'Check with S2 columns line');
   	is(scalar @header_string_s2_bits, scalar @exp2_bits,'Check with S2 correct number of lines');

   	foreach my $got_line (@header_string_s2_bits){
   		my $expected_line = shift @exp2_bits;
   		is_deeply(header_to_hash($got_line),header_to_hash($expected_line),'Checking with S2 header line');
   	}

  };

  subtest 'gen_record' => sub {

	my $contig_1 = new Sanger::CGP::Vcf::Contig(-name => '1', -length => 10, -species => 'jelly fish', -assembly => '25');
  	my $contig_x = new Sanger::CGP::Vcf::Contig(-name => 'x', -length => 10, -species => 'jelly fish', -assembly => '25');

  	my $contigs = [$contig_1,$contig_x];

  	my $converter = new Sanger::CGP::Pindel::OutputGen::VcfConverter(-contigs => $contigs);

	my $record_1 = new Sanger::CGP::Pindel::OutputGen::CombinedRecord(
	  -id => 'id1',
	  -idx => 'D12',
	  -reads => [],
	  -chro => 'x',
	  -start => 30,
	  -end => 35,
	  -range_start => 25,
	  -range_end => 40,
	  -length => 5,
	  -min_change => 'atg',
	  -lub => 'c',
	  -ref_seq => '',
	  -alt_seq => 'atg',
	  -sum_ms => 234,
	  -s1 => 3,
	  -s2 => 4,
	  -s2_presant => 1,
	  -type => 'I',
	  -repeats => 1,
	  -p_mt_pos => 1,
	  -p_mt_neg => 2,
	  -p_wt_pos => 3,
	  -p_wt_neg => 4,
	  -b_mt_pos => 5,
	  -b_mt_neg => 6,
	  -b_wt_pos => 7,
	  -b_wt_neg => 8,
	  -d_mt_pos => 9,
	  -d_mt_neg => 10,
	  -d_wt_pos => 11,
	  -d_wt_neg => 12,
	  -rd_mt_pos => 13,
	  -rd_mt_neg => 14,
	  -rd_wt_pos => 15,
	  -rd_wt_neg => 16,
	  -uc_mt_pos => 17,
	  -uc_mt_neg => 18,
	  -uc_wt_pos => 19,
	  -uc_wt_neg => 20,
	  -total_wt_rg_count => 5,
	  -call_wt_rg_count => 1,
	  -total_mt_rg_count => 9,
	  -call_mt_rg_count => 7,
	);

	my $exp1 = qq{x	30	id1	C	CATG	234	.	PC=I;RS=25;RE=40;LEN=5;S1=3;S2=4;REP=1	GT:PP:NP:PB:NB:PD:ND:PR:NR:PU:NU:TG:VG	./.:3:4:7:8:11:12:15:16:19:20:5:1	./.:1:2:5:6:9:10:13:14:17:18:9:7\n};

	is($converter->gen_record($record_1),$exp1,'gen_record');
  };
};


## breaks up a header line into manageable chunks for comparison. This way order shouldnt matter....
sub header_to_hash{
	my $line = shift;

	my ($k,$data) = $line =~ /^##([^=]+)=<(.+)>$/;

	unless(defined $k) {
	  ($k,$data) = $line =~ /^##([^=]+)=(.+)$/;
	}

	if($k =~ m/^vcfProcessLog_/) {
	  $k =~ s/_.*//;
	}

	my $ret = {key => $k};
	my $fragments = [];

	my @bits = split(/,/,$data);
	foreach my $pair (@bits){
		my($key,$val) = split /=</, $pair;
		if(defined $val) {
		  chop $val;
		}
		else {
		  ($key,$val) = split /=/, $pair;
		}

		unless($val){
			push @$fragments, $key;
		}else{
			$ret->{$key} = $val;
		}
	}

	$ret->{'frag'} = $fragments;
  return $ret;
}


done_testing();
