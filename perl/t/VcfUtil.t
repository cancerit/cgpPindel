use strict;
use Test::More;
use Test::Fatal;
use Data::Dumper;

my $MODULE = 'Sanger::CGP::Vcf::VcfUtil';

use DateTime qw();

use_ok 'Sanger::CGP::Vcf::Sample';
use_ok 'Sanger::CGP::Vcf::Contig';
use_ok 'Sanger::CGP::Vcf::VcfProcessLog';
use_ok 'Sanger::CGP::Vcf::VcfUtil';
use_ok 'Vcf';

subtest 'Non-object funcions' => sub {


  subtest 'add_vcf_sample' => sub {

  	my $wt_sample = new Sanger::CGP::Vcf::Sample(-name => 'test_wt');

    my $exp7 = q{##fileformat=VCFv4.1
##SAMPLE=<ID=test_wt,SampleName=test_wt>
##SAMPLE=<ID=test_wt_2,SampleName=test_wt_2,Study=test_study>
##SAMPLE=<ID=test_wt_3,SampleName=test_wt_3,Study=test_study,Protocol=protocol1>
##SAMPLE=<ID=test_wt_4,SampleName=test_wt_4,Study=test_study,Protocol=protocol1,Source=COSMIC>
##SAMPLE=<ID=test_wt_5,SampleName=test_wt_5,Description="Blah Blah Blah",Study=test_study,Source=COSMIC,Protocol=protocol1>
##SAMPLE=<ID=test_wt_6,SampleName=test_wt_6,Description="Blah Blah Blah",Study=test_study,Source=COSMIC,Accession=COS_1231,Protocol=protocol1>
##SAMPLE=<ID=test_wt_7,SampleName=test_wt_7,Description="Blah Blah Blah",Study=test_study,Source=COSMIC,Accession=COS_1231,Platform=plat1,Protocol=protocol1>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	test_wt	test_wt_2	test_wt_3	test_wt_4	test_wt_5	test_wt_6	test_wt_7
};

	my $vcf = new Vcf;

  	Sanger::CGP::Vcf::VcfUtil::add_vcf_sample($vcf,$wt_sample);

    $wt_sample->name('test_wt_2');
    $wt_sample->study('test_study');
  	Sanger::CGP::Vcf::VcfUtil::add_vcf_sample($vcf,$wt_sample);

    $wt_sample->name('test_wt_3');
  	$wt_sample->seq_protocol('protocol1');
  	Sanger::CGP::Vcf::VcfUtil::add_vcf_sample($vcf,$wt_sample);

  	$wt_sample->name('test_wt_4');
  	$wt_sample->accession_source('COSMIC');
  	Sanger::CGP::Vcf::VcfUtil::add_vcf_sample($vcf,$wt_sample);

  	$wt_sample->name('test_wt_5');
  	$wt_sample->description('Blah Blah Blah');
  	Sanger::CGP::Vcf::VcfUtil::add_vcf_sample($vcf,$wt_sample);

  	$wt_sample->name('test_wt_6');
  	$wt_sample->accession('COS_123');
  	Sanger::CGP::Vcf::VcfUtil::add_vcf_sample($vcf,$wt_sample);

  	$wt_sample->name('test_wt_7');
  	$wt_sample->platform('plat1');
  	Sanger::CGP::Vcf::VcfUtil::add_vcf_sample($vcf,$wt_sample);


  	my $header_samples = $vcf->format_header;
    my @header_samples_bits = split "\n", $header_samples;
   	my $header_samples_columns = pop @header_samples_bits;
   	my @exp1_bits = split "\n", $exp7;
   	my $exp7_columns = pop @exp1_bits;

   	is($header_samples_columns,$exp7_columns,'Check samples columns line');
   	is(scalar @header_samples_bits, scalar @exp1_bits,'Check correct number of lines');

   	foreach my $got_line (@header_samples_bits){
   		my $expected_line = shift @exp1_bits;
   		is_deeply(header_to_hash($got_line),header_to_hash($got_line),'Checking header line');
   	}
  	#is_deeply([split "\n", $vcf->format_header],[split "\n",$exp7],'Check Sample on header string');

  };



  subtest 'gen_tn_vcf_header' => sub {

  	my $wt_sample = new Sanger::CGP::Vcf::Sample(-name => 'test_wt');
  	my $mt_sample = new Sanger::CGP::Vcf::Sample(-name => 'test_mt');

  	my $contig_1 = new Sanger::CGP::Vcf::Contig(-name => '1', -length => 10, -species => 'jelly fish', -assembly => '25');
  	my $contig_x = new Sanger::CGP::Vcf::Contig(-name => 'x', -length => 10, -species => 'jelly fish', -assembly => '25');

  	my $contigs = [$contig_1,$contig_x];

  	my @process_logs = ();

    push @process_logs, new Sanger::CGP::Vcf::VcfProcessLog(
	  -input_vcf_source => 'banana',
      -input_vcf_ver => 'test_version1',
	);

	push @process_logs, new Sanger::CGP::Vcf::VcfProcessLog(
	  -input_vcf_source => 'Jeff K',
      -input_vcf_ver => 'test_version2',
      -input_vcf_params => {'o'=>'my_file'},
	);

  	my $reference_name = 'BOB_ref';
    my $input_source = 'BOB';
  	my $date = DateTime->now->strftime('%Y%m%d');

	my $format = [
		{key => 'FORMAT', ID => 'GT', Number => 1, Type => 'String', Description => 'Genotype'},
		{key => 'FORMAT', ID => 'PP', Number => 1, Type => 'Integer', Description => 'Pindel calls on the positive strand'},
	];

	my $info = [
		{key => 'INFO', ID => 'PC', Number => 1, Type => 'String', Description => 'Pindel call'},
		{key => 'INFO', ID => 'RS', Number => 1, Type => 'Integer', Description => 'Range start'},
	];

	my $other = [
		{key => 'random', value => 'This is random'},
		{key => 'VERY_RANDOM', value => 'So is this'},
	];

  	my $exp1 =
qq{##fileformat=VCFv4.1
##fileDate=$date
##source_$date.1=BOB
##reference=BOB_ref
##contig=<ID=1,assembly=25,length=10,species="jelly fish">
##contig=<ID=x,assembly=25,length=10,species="jelly fish">
##INFO=<ID=PC,Number=1,Type=String,Description="Pindel call">
##INFO=<ID=RS,Number=1,Type=Integer,Description="Range start">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PP,Number=1,Type=Integer,Description="Pindel calls on the positive strand">
##random=This is random
##VERY_RANDOM=So is this
##vcfProcessLog=<InputVCFSource=<banana>,InputVCFVer=<test_version1>>
##vcfProcessLog=<InputVCFSource=<Jeff K>,InputVCFVer=<test_version2>,InputVCFParam=<o=my_file>>
##SAMPLE=<ID=NORMAL,SampleName=test_wt>
##SAMPLE=<ID=TUMOUR,SampleName=test_mt>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NORMAL	TUMOUR
};
	is_deeply([split "\n", Sanger::CGP::Vcf::VcfUtil::gen_tn_vcf_header($wt_sample, $mt_sample, $contigs, \@process_logs, $reference_name, $input_source, $info, $format, $other, "\n")],[split "\n",$exp1],'Header string construction');

  };

};

## breaks up a header line into manageable chunks for comparison. This way order shouldnt matter....
sub header_to_hash{
	my($line) = @_;

	my ($k,$data) = $line =~ /^##(.+)=<()>^/;

	my $ret = {key => $k};


	my @bits = split(/,/,$data);

	my $fragments = [];

	my $key_gen = 0;

	foreach my $pair (@bits){
		my($key,$val) = split /=/, $pair;

		unless($val){
			push @$fragments, $fragments;
		}else{
			$ret->{$key} = $val;
		}
	}

	$ret->{'frag'} = $fragments;
}

done_testing();
