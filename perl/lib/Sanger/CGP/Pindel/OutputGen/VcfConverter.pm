package Sanger::CGP::Pindel::OutputGen::VcfConverter;

########## LICENCE ##########
# Copyright (c) 2014-2017 Genome Research Ltd.
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


use Sanger::CGP::Pindel;

use strict;
use Sanger::CGP::Vcf::VcfUtil;

use constant SEP => "\t";
use constant NL => "\n";

1;

sub new{
	my $proto = shift;
	my (%args) = @_;
	my $class = ref($proto) || $proto;

	my $self = {};
	bless $self, $class;

	$self->init(%args);

	return $self;
}

sub init{
	my($self,%args) = @_;
	$self->{_contigs} = $args{-contigs},
	$self->{_format} = 'GT:PP:NP:PB:NB:PD:ND:PR:NR:PU:NU';
}


=head init_tn_vcf_header

Generates a Vcf header String for NORMAL/TUMOUR comparisons.

@param1 wt_sample      - a Sanger::CGP::Vcf::Sample object representing the wild type sample.

@param2 mt_sample      - a Sanger::CGP::Vcf::Sample object representing the mutant type sample.

@param3 process_logs   - an array-ref of Sanger::CGP::Vcf::VcfProcessLog objects.

@param4 s2_presant     - boolean. Set to 1 if the S2 column is presant in the Vcf input file(s).

@param5 reference_name - a String containing the name of the reference used in the VCF.

@param6 input_source   - a String containing the name and version of the application or source of the VCF data.

=cut
sub gen_header{
	my($self,$wt_sample, $mt_sample, $process_logs, $s2_presant, $reference_name, $input_source) = @_;

	my $contigs = $self->{_contigs};

	my $info = [
		{key => 'INFO', ID => 'PC', Number => 1, Type => 'String', Description => 'Pindel call'},
		{key => 'INFO', ID => 'RS', Number => 1, Type => 'Integer', Description => 'Range start'},
		{key => 'INFO', ID => 'RE', Number => 1, Type => 'Integer', Description => 'Range end'},
		{key => 'INFO', ID => 'LEN', Number => 1, Type => 'Integer', Description => 'Length'},
		{key => 'INFO', ID => 'REP', Number => 1, Type => 'Integer', Description => 'Change repeat count within range'},
		{key => 'INFO', ID => 'S1', Number => 1, Type => 'Integer', Description => 'S1'}
	];

	push @{$info}, {key => 'INFO', ID => 'S2', Number => 1, Type => 'Float', Description => 'S2'} if($s2_presant);

	my $format = [
		{key => 'FORMAT', ID => 'GT', Number => 1, Type => 'String', Description => 'Genotype'},
		{key => 'FORMAT', ID => 'PP', Number => 1, Type => 'Integer', Description => 'Pindel calls on the positive strand'},
		{key => 'FORMAT', ID => 'NP', Number => 1, Type => 'Integer', Description => 'Pindel calls on the negative strand'},
		{key => 'FORMAT', ID => 'PB', Number => 1, Type => 'Integer', Description => 'BWA calls on the positive strand'},
		{key => 'FORMAT', ID => 'NB', Number => 1, Type => 'Integer', Description => 'BWA calls on the negative strand'},
		{key => 'FORMAT', ID => 'PD', Number => 1, Type => 'Integer', Description => 'BWA mapped reads on the positive strand'},
		{key => 'FORMAT', ID => 'ND', Number => 1, Type => 'Integer', Description => 'BWA mapped reads on the negative strand'},
		{key => 'FORMAT', ID => 'PR', Number => 1, Type => 'Integer', Description => 'Total mapped reads on the positive strand'},
		{key => 'FORMAT', ID => 'NR', Number => 1, Type => 'Integer', Description => 'Total mapped reads on the negative strand'},
		{key => 'FORMAT', ID => 'PU', Number => 1, Type => 'Integer', Description => 'Unique calls on the positive strand'},
		{key => 'FORMAT', ID => 'NU', Number => 1, Type => 'Integer', Description => 'Unique calls on the negative strand'},
	];

	return Sanger::CGP::Vcf::VcfUtil::gen_tn_vcf_header( $wt_sample, $mt_sample, $contigs, $process_logs, $reference_name, $input_source, $info, $format, []);
}

sub gen_record{
	my($self, $record) = @_;

	# CHR POS ID REF ALT QUAL FILTER INFO FORMAT GENO GENO

	my $start = $record->start();
	$start-- if(substr($record->type(),0,1) eq 'D');

	my $ret = $record->chro().SEP;
	$ret .= $start.SEP;
	$ret .= $record->id().SEP;

	my $ref = uc ($record->lub . $record->ref_seq);
	my $alt = uc ($record->lub . $record->alt_seq);

	$ret .= $ref.SEP;
	$ret .= $alt.SEP;
	$ret .= $record->sum_ms().SEP;
	$ret .= '.'.SEP;

	# INFO
	#PC=D;RS=19432;RE=19439;LEN=3;S1=4;S2=161.407;REP=2;PRV=1
	$ret .= 'PC='.$record->type().';';
	$ret .= 'RS='.$record->range_start().';';
	$ret .= 'RE='.$record->range_end().';';
	$ret .= 'LEN='.$record->length().';';
	$ret .= 'S1='.$record->s1().';';
	$ret .= 'S2='.$record->s2().';' if(defined $record->s2()); ## not presant in older versions of pindel
	$ret .= 'REP='.$record->repeats().SEP;

	# FORMAT
	$ret .= $self->{_format}.SEP;

	# 'GT:PP:NP:PB:NB:PD:ND:PR:NR:PU:NU'
	# NORMAL GENO
	$ret .= './.:';
	$ret .= $record->p_wt_pos().':';
	$ret .= $record->p_wt_neg().':';
	$ret .= $record->b_wt_pos().':';
	$ret .= $record->b_wt_neg().':';
	$ret .= $record->d_wt_pos().':';
	$ret .= $record->d_wt_neg().':';
	$ret .= $record->rd_wt_pos().':';
	$ret .= $record->rd_wt_neg().':';
	$ret .= $record->uc_wt_pos().':';
	$ret .= $record->uc_wt_neg().SEP;

	# TUMOUR GENO
	$ret .= './.:';
	$ret .= $record->p_mt_pos().':';
	$ret .= $record->p_mt_neg().':';
	$ret .= $record->b_mt_pos().':';
	$ret .= $record->b_mt_neg().':';
	$ret .= $record->d_mt_pos().':';
	$ret .= $record->d_mt_neg().':';
	$ret .= $record->rd_mt_pos().':';
	$ret .= $record->rd_mt_neg().':';
	$ret .= $record->uc_mt_pos().':';
	$ret .= $record->uc_mt_neg().NL;

	return $ret;
}
