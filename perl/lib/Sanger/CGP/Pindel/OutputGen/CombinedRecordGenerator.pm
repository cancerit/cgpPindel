package Sanger::CGP::Pindel::OutputGen::CombinedRecordGenerator;

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


use Sanger::CGP::Pindel;

use strict;
use Carp;
use English qw( -no_match_vars );
use Data::Dumper;

use Bio::DB::Sam;

use Sanger::CGP::Pindel::OutputGen::CombinedRecord;

use base qw(Sanger::CGP::Pindel::OutputGen::PindelRecordParser);

use Const::Fast qw(const);

const my $F_UNMAPPED    => 4;
const my $F_NOT_PRIMARY => 256;
const my $F_QC_FAIL     => 512;
const my $F_DUPLICATE   => 1024;
const my $F_SUPP_ALIGN  => 2048;

#ordered most-frequent to least for efficiency
const my @FILTER_READS_IF_FLAG => ($F_DUPLICATE, $F_UNMAPPED, $F_SUPP_ALIGN, $F_QC_FAIL, $F_NOT_PRIMARY);

sub new{
	my ($class, %args) = @_;
	my $self = $class->SUPER::new(%args);
	return $self;
}

sub init{
	my($self,%args) = @_;
	croak "Undefined wt_sam obj argument" unless(defined $args{-wt_sam});
	croak "Undefined mt_sam obj argument" unless(defined $args{-mt_sam});
	$self->{_wt_sam} = $args{-wt_sam};
	$self->{_mt_sam} = $args{-mt_sam};
	$self->{_mutant_sample_name} = $args{-mutant_sample_name};
	$self->SUPER::init(%args);## this has to be at the end because we are overriding the next record method with our own
}

=head2 next_record

Grabs the next record from the pindel output file.
Overrides next_record in 'Sanger::CGP::Pindel::OutputGen::PindelRecordParser'.

@returns - a Sanger::CGP::Pindel::OutputGen::CombinedRecord;

=cut
sub next_record{
	my($self) = @_;
	my $ret = $self->{_current_record};
	$self->{_current_record} = $self->_process_record(new Sanger::CGP::Pindel::OutputGen::CombinedRecord());
	return $ret;
}

=head2 _process_record

Overrides _process_record in 'Sanger::CGP::Pindel::OutputGen::PindelRecordParser'.

@param1 record - requires an empty Sanger::CGP::Pindel::OutputGen::CombinedRecord object.

@returns       - the Sanger::CGP::Pindel::OutputGen::CombinedRecord

=cut
sub _process_record{
	my ($self,$record) = @_;

	# Process the pindel part of the record.
	return undef unless $self->SUPER::_process_record($record);

	$self->_process_counts($record,'mt');
	$self->_process_counts($record,'wt');

	return $record;
}

=head2 _process_counts

Use the internal sam object to generate counts for the BWA calls and depth for a given record.

@param1 record    - A Sanger::CGP::Pindel::OutputGen::CombinedRecord object.

@param2 samp_type - A string, either 'wt' or 'mt'.

@returns the original Sanger::CGP::Pindel::OutputGen::CombinedRecord.

=cut
sub _process_counts{
	my ($self, $record, $samp_type) = @_; #wt mt...'

	# Grab the correct sam object from self...
	my $sam_obj = $self->{"_${samp_type}_sam"};

	croak "Sam object not defined _${samp_type}_sam." unless $sam_obj;

	## BWA Depths... not taking into account reads passed to Pindel...
	my ($pos_bwa_reads, $neg_bwa_reads);
	if($record->type eq 'I'){
		($pos_bwa_reads, $neg_bwa_reads) = _read_depth($sam_obj, $record->chro, $record->range_start, $record->range_end-1); # as indels are always attached to the base before in Bio::DB::Sam
	}else{
		($pos_bwa_reads, $neg_bwa_reads) = _read_depth($sam_obj, $record->chro, $record->range_start, $record->range_start);
	}

	# b_wt_pos this is the method name structure of CombinedRecord methods.

	# bwa calls
	my($pos_call_reads, $neg_call_reads) = _count_sam_event_reads($record, $sam_obj, $samp_type);

	## set up the correct wt/mt method names for the record
	## TODO might have to revisit this bit - perhaps split the counts by sample
	my $b_pos_method = "b_${samp_type}_pos";
	my $b_neg_method = "b_${samp_type}_neg";
	my $d_pos_method = "d_${samp_type}_pos";
	my $d_neg_method = "d_${samp_type}_neg";
	my $p_pos_method = "p_${samp_type}_pos";
	my $p_neg_method = "p_${samp_type}_neg";
	my $rd_pos_method = "rd_${samp_type}_pos";
	my $rd_neg_method = "rd_${samp_type}_neg";
	my $uc_pos_method = "uc_${samp_type}_pos";
	my $uc_neg_method = "uc_${samp_type}_neg";
	my $total_rg_count_method = "total_${samp_type}_rg_count";
	my $call_rg_count_method = "call_${samp_type}_rg_count";

	# bwa depth
	$record->$d_pos_method(scalar keys %$pos_bwa_reads);
	$record->$d_neg_method(scalar keys %$neg_bwa_reads);

	## ok so... we have all the read names from the pass through pindel.
	## we also have all the read names from the bam file...
	## thus we should be able to get unique counts.....

	my %pos_pin_reads = ();
	my %neg_pin_reads = ();

	foreach my $pindel_sample ($record->samples()){
		my $pin_sample_type = $pindel_sample eq $self->{_mutant_sample_name} ? 'mt' : 'wt';

		## this method is already tied to a sample type we need to skip mismatching sampes from the list. As we do not store the normal we still need to loop.
		next unless $pin_sample_type eq $samp_type;

		my $pos_reads = $record->get_reads($pindel_sample,'+') || [];
		my $neg_reads = $record->get_reads($pindel_sample,'-') || [];

		## In the pindel parser we put a chunk of text on the end of the read to make it unique in cases where it might be used in more than one variant.
		## This needs to be removed so that we can merge on unique read name....
		## TODO the read groups might not be needed....
		foreach my $read (@$pos_reads){

			my ($rg) = grep {/RG:Z:.+/} @$read;
			$rg = '' unless $rg;
			$rg =~ s/^RG:Z://;

			my $name = $read->[0];
			$name =~ s/_[^_]+_[^_]+$//;
			$pos_pin_reads{$name} = $rg;
		}

		foreach my $read (@$neg_reads){

			my ($rg) = grep {/RG:Z:.+/} @$read;
			$rg = '' unless $rg;
			$rg =~ s/^RG:Z://;

			my $name = $read->[0];
			$name =~ s/_[^_]+_[^_]+$//;
			$neg_pin_reads{$name} = $rg;
		}

		## merge the hashes together... dont worry about the values as they are not important.
		## This is an efficient way of hash merging
		## Depth first..

		@{$pos_bwa_reads}{keys %pos_pin_reads} = values %pos_pin_reads;
		@{$neg_bwa_reads}{keys %neg_pin_reads} = values %neg_pin_reads;

		## Then calls...
		@{$pos_call_reads}{keys %pos_pin_reads} = values %pos_pin_reads;
		@{$neg_call_reads}{keys %neg_pin_reads} = values %neg_pin_reads;
	}

	## Read group counts
	my %unique_rg;
	$unique_rg{$_}++ for (values %$pos_bwa_reads,values %$neg_bwa_reads);

	my %unique_call_rg;
	$unique_call_rg{$_}++ for (values %$pos_call_reads,values %$neg_call_reads);

	# read group counts
	$record->$total_rg_count_method(scalar keys %unique_rg);
	$record->$call_rg_count_method(scalar keys %unique_call_rg);

	# pindel calls
	$record->$p_pos_method(scalar keys %pos_pin_reads);
	$record->$p_neg_method(scalar keys %neg_pin_reads);

	# real depth
	$record->$rd_pos_method(scalar keys %$pos_bwa_reads);
	$record->$rd_neg_method(scalar keys %$neg_bwa_reads);

	# unique calls
	$record->$uc_pos_method(scalar keys %$pos_call_reads);
	$record->$uc_neg_method(scalar keys %$neg_call_reads);

	return $record;
}

=head2 _read_depth

Performs a pileup for a given location and returns a hash of al the reads and their read groups by strand.

@param1 sam    - A Bio::DB:Sam object.

@param2 seqid  - A string representing a chromosome or contig used in the reference.

@param3 start  - An integer genomic start coordinate.

@param4 end    - An integer genomic end coordinate.

@returns (Hash-ref->{read_name}=>read_group, Hash_ref->{read_name}=>read_group) positive reads, negative reads.

=cut
sub _read_depth {
	my ($sam, $seqid, $r_start, $r_end) = @_;

	my $align_iter = $sam->get_features_by_location(
		-seq_id => $seqid,
		-start => $r_start,
		-end => $r_end,
		-filter => \&_read_filter,
		-iterator => 1,
	);
	my %pos_loc_tmp;
	my %neg_loc_tmp;

  while(my $a = $align_iter->next_seq) {
#	foreach my $a(@alignments) {
		my $rn = $a->qname;

		my $rg = $a->aux_get('RG');
		$rg = '' unless $rg;

		$rn =~ s/_r[0-9]+_[DI]{1,2}[0-9]+$//; # where the read has come from a pindel bam file.... TODO this might not be needed any more...
		if($a->reversed) {
			$neg_loc_tmp{$rn} = $rg;
		}else {
			$pos_loc_tmp{$rn} = $rg;
		}
	}

	return (\%pos_loc_tmp,\%neg_loc_tmp);
}

sub _read_filter {
  my $flag = shift->flag;
  for my $test(@FILTER_READS_IF_FLAG) {
    return 0 if(($flag | $test) == $flag);
  }
  return 1;
}

sub _count_sam_event_reads{
	my ($record, $sam_obj, $samp_type_key) = @_; #wt mt...'

	## counts [0] = negative strand [1] = positive strand
	my @ins_count = (0,0);
	my @del_count = (0,0);

	# an array of hashes here we want to collect all the read names
	# reads [0] = negative strand [1] = positive strand
	my @ins_reads = ({},{});
	my @del_reads = ({},{});

	## BWA event counts..... but why not DI????? There is no representation of DI in BWA pileups.... counts were added with UC correction...
	## Unless the event is a DI collect all of the events spotted within the range start/end

	# These will be used by the pileup callback.
	my $g_r_start = $record->range_start;
	my $g_r_end = $record->range_end-1; # as indels are always attached to the base before in Bio::DB::Sam
	my $g_chg_len = $record->type eq 'D' ? length $record->ref_seq : length $record->alt_seq;
	my $g_min_chg_len = length $record->min_change;
	my $g_e_type = $record->type;
	my $chr = $record->chro;
	my $g_match = '';

	## DIs dont really exist as BWA calls more likely to be del followed by an ins...
	if($g_e_type eq 'DI') {
		my $g_i_len = $g_chg_len;
		$g_chg_len = $record->end - $record->start; # length of the deletion
		$g_match = $g_chg_len.'D'.$g_i_len.'I';
	}

    croak 'Change length is less than 1' if $g_chg_len == 0;

    ## define all the variables outside of the loop to save on memory allocation events..
    my ($seqid,$pos,$pileup); ##1
    my ($pu, $value, $abs, $strand_index, $rn); ##2

	$sam_obj->fast_pileup("$chr:$g_r_start-$g_r_end", sub {
			($seqid,$pos,$pileup) = @_; ##1

			if($pos >= $g_r_start && $pos <= $g_r_end) {
				##2
				## no critic
				foreach $pu (@{$pileup}) {
 					$value = $pu->indel;
        ## use critic
 					next if($value == 0); ## we are not interested in non-call reads...
 					next if($g_e_type eq 'D'  && $value > 0);
					next if($g_e_type eq 'I'  && $value < 0);

					my $a = $pu->alignment;
					next unless(_read_filter($a));

					next if($g_e_type eq 'DI' && $a->cigar_str !~ $g_match);
					$abs = abs $value;

				#	warn "change len: $abs $g_chg_len, ".$record->ref_seq;

					##TODO WE SHOULD LOOK AT MATCHING ON CHANGE rather than length????????????? jwh
					if($abs == $g_chg_len || $abs == $g_min_chg_len) {
						$strand_index = $a->reversed == 1 ? 0 : 1;
 					    $rn = $a->qname;
						$rn =~ s/_r[0-9]+_[DI]{1,2}[0-9]+$//;

						if($value > 0) {
							$ins_count[$strand_index]++;
							$ins_reads[$strand_index]->{$rn}++;
						}else{
							$del_count[$strand_index]++;
							$del_reads[$strand_index]->{$rn}++;
						}
					}
				}
			}
			return;
    	}
	);

	my $b_pos_method = "b_${samp_type_key}_pos";
	my $b_neg_method = "b_${samp_type_key}_neg";

	## Set the bwa calls
	# b_wt_pos this is the method name structure of CombinedRecord
	if($g_e_type eq 'D'){
		$record->$b_pos_method($del_count[1]);
		$record->$b_neg_method($del_count[0]);
	}elsif($g_e_type eq 'I'){
		$record->$b_pos_method($ins_count[1]);
		$record->$b_neg_method($ins_count[0]);
	}elsif($g_e_type eq 'DI'){
		##TODO if the cigar matches the cigar pattern will have counts (this bit was added as part of the uc correction) not sure why zero... jwh
		$record->$b_pos_method(0);
		$record->$b_neg_method(0);
	}

	## Reorganise the reads into pos/neg strand groups. Outside of this method we do not care about ins/dels counts...
	my %pos_reads = %{$ins_reads[1]};
	my %neg_reads = %{$ins_reads[0]};

	@pos_reads{keys %{$del_reads[1]}} = values{%{$del_reads[1]}};
	@neg_reads{keys %{$del_reads[0]}} = values{%{$del_reads[0]}};

	## return the reads to the outside world, we need to merge these with the reads discovered from Pindel to get accurate counts...
	return (\%pos_reads,\%neg_reads);
}

1;
