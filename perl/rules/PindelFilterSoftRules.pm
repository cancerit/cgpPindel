####################################################
# Copyright (c) 2013 Genome Research Ltd.
# Author: Cancer Genome Project, cgpit@sanger.ac.uk
# See LICENCE.TXT for details
####################################################
{
	tag  => 'INFO/LEN',                   # The VCF tag to apply this filter on
	name => 'F017',                       # The filter ID
	desc => 'Variant must not overlap with a simple repeat',  # Description for the VCF header
	test => sub {

#update PINDEL_VARIANT pvo
#set pvo.FLAG_RUN_BIT_CODE = bitor(pvo.FLAG_RUN_BIT_CODE,?)
#,pvo.flag_bit_code = nvl(
#(SELECT bitand(pvo.FLAG_BIT_CODE, bitnot(?))
#FROM REFERENCE_REPEAT rr
#WHERE rr.ID_SEQUENCE = pvo.ID_SEQUENCE
#AND
#(
#	(rr.UNIT_LENGTH <= 2 AND rr.UNIT_LENGTH <= 10 AND rr.REPETITIONS >= 8)
#	OR
#	(rr.UNIT_LENGTH > 2 AND rr.UNIT_LENGTH <= 10 AND rr.REPETITIONS >= 5)
#)
#AND rr.REPEAT_TYPE = ?
#AND rr.MIN_POS <= (pvo.range_start + length(pvo.change))
#AND rr.MAX_POS >= (pvo.range_end - length(pvo.change))
#group by rr.ID_SEQUENCE
#),bitor(PVO.flag_bit_code,?))
#where pvo.id_pindel_run = ?
#and bitand(pvo.FLAG_BIT_CODE, ?) = ?
		#22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	1:0:1:0:4:5:4:6:1:1	3:1:4:2:25:20:32:21:18:8

		### HACK Dirty dirty dirty......
		unless($main::VCF_FLAGGING_REPEATS_TABIX){
			use Tabix;
			$main::VCF_FLAGGING_REPEATS_TABIX = new Tabix(-data => $ENV{VCF_FLAGGING_REPEATS},-index => $ENV{VCF_FLAGGING_REPEATS}.'.tbi');
		}

		my($from) = ";$$RECORD[7]" =~ m/;RS=(\d+)/;
		my($to) = ";$$RECORD[7]" =~ m/;RE=(\d+)/;

		my $ret = eval{
			my $res = $main::VCF_FLAGGING_REPEATS_TABIX->query($CHROM,($from-1),$to);
			return $PASS if(!defined $res->get); # no valid entries (chromosome not in index) so must pass
			return $FAIL if($main::VCF_FLAGGING_REPEATS_TABIX->read($res));
			return $PASS;
		};
		if($@) {
	    die $@;
		}
		return $ret;
	}
},
