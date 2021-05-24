package Sanger::CGP::Pindel::OutputGen::PindelRecordParser;

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


use Sanger::CGP::Pindel;

use strict;
use Carp;
use English qw( -no_match_vars );
use Data::Dumper;
use Const::Fast qw(const);

use Sanger::CGP::Pindel::OutputGen::PindelRecord;

const my $MAX_REPEAT_UNIT_SIZE => 5;

1;

sub new{
	my ($class, %args) = @_;
	my $self = { };
	bless $self, $class;
	$self->init(%args);
	return $self;
}

sub init{
	my($self,%args) = @_;

	my $fh = $args{-fh};

	unless (defined $fh){
		croak "Undefined file handle or file path argument" unless(defined $args{-path});
		open($fh,"<", $args{-path}) or croak "Cannot open |".$args{-path}."| for reading: $!";
		$self->{_path} = $args{-path};
	}

	$self->{_fh} = $fh;
	$self->{_fai} = $args{-fai};
	$self->{_noreads} = $args{-noreads} || 0;

	## clear the first line of ##+....

  my $first_line = <$fh>;
	if(defined $first_line){
		chomp $first_line;
		croak "Expecting the first line to be ##+.... |$first_line|" unless $first_line =~ m/^##+/;
	}

	$self->next_record();
}

sub close{
	my ($self) = @_;

	if($self->{_path} && defined $self->{_fh}){
		close $self->{_fh} or croak "Unable to close |".$self->{_path}."|";
	}

	$self->{_fh} = undef;
	$self->{_path} = undef;
	$self->{_version} = undef;
}

sub version{
	return shift->{_version};
}

=head3 next_record

Grabs the next record from the pindel output file.

@returns - a Sanger::CGP::Pindel::OutputGen::PindelRecord;

=cut
sub next_record{
	my($self) = @_;
	my $ret = $self->{_current_record};
	$self->{_current_record} = $self->_process_record(new Sanger::CGP::Pindel::OutputGen::PindelRecord());
	return $ret;
}

sub _process_record{
	my ($self,$record) = @_;

	my $fh = $self->{_fh};
	my $record_header = <$fh>;
	return undef unless(defined $record_header);
	chomp $record_header;

	my $ref_line = <$fh>;
	return undef unless(defined $ref_line);
	chomp $ref_line;

	## There is the old version and the new version of the output but there is no imediate way to tell
	#v1	Supports 3	+ 0	- 3	S1
	#v2	Supports 3      2       + 0     0       - 3     2 S1
	if($record_header =~ m/Supports\s\d+\s[+-]\s\d+\s[+-]\s\d+\sS1\s/){
		print STDERR 'File contains calls of different pindel versions\n' if(defined $self->{_version} && $self->{_version}  ne 'v01');
		$self->{_version} = 'v01';
		$record->version($self->{_version});
		_parse_header($record,$record_header);
	}else{
		print STDERR 'File contains calls of different pindel versions\n' if(defined $self->{_version} && $self->{_version}  ne 'v02');
		$self->{_version} = 'v02';
		$record->version($self->{_version});
		_parse_header_v02($record,$record_header);
	}



	my $alignments = [];

	# Collect all the reads...
	while(my $line = <$fh>) {
		last if $line =~  m/^##+/;
		chomp $line;

		push(@$alignments, $line);
	}

	$self->_parse_alignment($record, $alignments, \$ref_line);

	# this is where we can add additional information about GC content
	my $chr = $record->chro;
	my $lhs_end = $record->range_start;
	my $rhs_start = $record->range_end;
	my $fai = $self->{_fai};
	my $r_type = $record->type;

	my $seq = $fai->fetch(sprintf '%s:%d-%d', $chr, $lhs_end - 199, $lhs_end);
	my $lhs_gc = ($seq =~ tr/GCgc//)/200;

	$seq = $fai->fetch(sprintf '%s:%d-%d', $chr, $rhs_start, $rhs_start + 199);
	my $rhs_gc = ($seq =~ tr/GCgc//)/200;
	my $ref_tmp;
	if($record->type eq 'D'){
		$ref_tmp = $record->ref_seq;
	}
	else {
		$ref_tmp = $record->alt_seq;
	}
	my $rng_gc = ($ref_tmp =~ tr/GCgc//) / length $ref_tmp;
	$record->gc_5p($lhs_gc);
	$record->gc_3p($rhs_gc);
	$record->gc_rng($rng_gc);

	return $record;
}

=head3 _parse_header

Processes the header component of a pindel record.

@param1 record         - a empty Sanger::CGP::Pindel::OutputGen::PindelRecord object.

@param2 alignment      - a string containing the pindel call header line.

=cut
sub _parse_header {
	my ($record, $head_line) = @_;

	my @head_bits = split /\t/, $head_line;
	my $call_num = shift @head_bits;

	# sort out event types
	my ($type, $event_size) = split / /, shift @head_bits;

	my $nt_field = shift @head_bits; # gives the non-templated sequence
	my ($nt_size,$nt_seq) = $nt_field =~ m/NT (\d+) "(.*)"/;

	## mark the variant as complex...
	if($type eq 'D' && $nt_size != 0) {
		$type .= 'I';
	}

	$record->length($event_size);
	$record->type($type);
	$record->idx($type.$call_num);
	$record->alt_seq($nt_seq);

	# set up coordinates
	my (undef, $chro) = split / /, shift @head_bits;
	my (undef, $start) = split / /, shift @head_bits;
	my $end = shift @head_bits;

	# moved this up so that DI range start and end in the same maner as D and I... jwh
	if($type eq 'DI') {
		$record->range_start($start);
		$record->range_end($end);
	}else {
		$record->range_start((split / /, shift @head_bits)[1]);
		$record->range_end(shift @head_bits);
	}

	## Pindel reports the bases either side of the event. This is natural for insertions but a bit obscure for deletions.
	$start++ unless ($type eq 'I');
	$end--; ## This brings the end of deletions back to the seq that is being deleted. This also makes insertions a single coordinate

	$record->chro($chro);
	$record->start($start);
	$record->end($end);

	# ignore these summaries, get direct from alignments
	shift @head_bits; # supporting reads
	shift @head_bits; # +ve supporting reads
	shift @head_bits; # -ve suporting reads

	# populate scores
	$record->s1((split / /, shift @head_bits)[1]);
	$record->s2((split / /, shift @head_bits)[1]) if($type ne 'DI') ;
	$record->sum_ms((split / /, shift @head_bits)[1]);
	$record->num_samples((split / /, shift @head_bits)[1]);
	my %sample_contrib;
	while(scalar @head_bits) {
	  my ($sample, $count) = split / /, shift @head_bits;
	  $sample_contrib{$sample} = $count;
	}
	$record->sample_contrib(\%sample_contrib);

	# remainder of data is obtained from aligments
	return $record;
}

=head3 _parse_header_v02

Processes the header component of a pindel record.

@param1 record    - a empty Sanger::CGP::Pindel::OutputGen::PindelRecord object.

@param2 alignment - a string containing the pindel call header line.

=cut
sub _parse_header_v02 {
	my ($record, $head_line) = @_;

	my @head_bits = split /\t/, $head_line;
	my $call_num = shift @head_bits;

	# sort out event types
	my ($type, $event_size) = split / /, shift @head_bits;

	my $nt_field = shift @head_bits; # gives the non-templated sequence
	my ($nt_size,$nt_seq) = $nt_field =~ m/NT (\d+) "(.*)"/;

	if($type eq 'D' && $nt_size != 0) {
		$type .= 'I';
	}

	$record->length($event_size);
	$record->type($type);
	$record->idx($type.$call_num);
	$record->alt_seq($nt_seq);

	# set up coordinates
	my (undef, $chro) = split / /, shift @head_bits;
	my (undef, $start) = split / /, shift @head_bits;
	my $end = shift @head_bits;

	## Pindel reports the bases either side of the event. This is natural for insertions but a bit obscure for deletions.
	$start++ unless ($type eq 'I');
	$end--; ## This brings the end of deletions back to the seq that is being deleted. This also makes insertions a single coordinate

	$record->chro($chro);
	$record->start($start);
	$record->end($end);

	if($type eq 'DI') {
		$record->range_start($start);
		$record->range_end($end);
	}else {
		$record->range_start((split / /, shift @head_bits)[1]);
		$record->range_end(shift @head_bits);
	}

	# ignore these summaries, get direct from alignments
	shift @head_bits; # supporting reads
	shift @head_bits; # supporting unique reads
	shift @head_bits; # +ve supporting reads
	shift @head_bits; # +ve supporting unique reads
	shift @head_bits; # -ve suporting reads
	shift @head_bits; # -ve suporting unique reads

	# populate scores
	$record->s1((split / /, shift @head_bits)[1]);
	$record->sum_ms((split / /, shift @head_bits)[1]);

	# remainder of data is obtained from aligments
	return $record;
}


=head3 _parse_alignment

processes the alignment component of a pindel record.

@param1 record         - a Sanger::CGP::Pindel::OutputGen::PindelRecord object. This must have the following fields
                         filled in: type, idx, chro, start, ref_seq, alt_seq.

@param2 alignment      - an array-ref of pindel alignment strings.

@param3 ref_line       - a String_ref of the pindel reference string.

=cut
sub _parse_alignment {
	my ($self, $record, $alignment, $ref_line) = @_;

	my ($ref_left, $ref_change, $ref_right) = ($$ref_line =~ m/([A-Z]+)(\s+|[a-z]+|[a-z]+.*[a-z]+)([A-Z]+)/);

	$record->ref_left($ref_left);
	$record->ref_right($ref_right);

	my $change_ref_offset = length $ref_left;
	my $change_ref_offset_end = $change_ref_offset + length $ref_change;
	my $record_type = $record->type();
    my $chr = $record->chro;
    my $record_idx = '_'.$record->idx;
    my $fai = $self->{_fai};

	# In the case of large deletions pindel replaces the deleted seq with a range <240> because of this we need to collect
	# the deleted seq from the reference. Also DIs do not display the deleted seq.....
	if($record_type eq 'I') {
	  $ref_change = q{}; ## $ref_change =~ s/\s//; ##remove all white space as there is no deleted ref_seq
	}elsif($record_type eq 'DI' || $ref_change =~ m/[[:digit:]]+/) {
		my $start = $record->start || $record->range_start;
		my $end = $record->end || $record->range_end;
		$ref_change = $fai->fetch("$chr:$start-$end");
	}

	$record->ref_seq($ref_change);
	$record->lub(substr($ref_left,-1)); # last unmodified base (useful for VCF conversion).

	# This is used to work out the the number of repeats also used alswhere.....
	if($record->alt_seq){
		$record->min_change(_shrink_change($record->alt_seq,$MAX_REPEAT_UNIT_SIZE)); # Shrink the change down to its minimum repeat component.
	}else{
		$record->min_change(_shrink_change($ref_change,$MAX_REPEAT_UNIT_SIZE)); # Shrink the change down to its minimum repeat component.
	}

	my $read_num = 1;
	my $ref_seq_length = length $ref_change;
	my $start_pos = $record->start();
	$start_pos++ if($record_type eq 'I');

	## Se up a local buffered region of reference. Rather than hitting the .fa file we can simply substr chunks out of this.....
	$self->{_buffer_region_chr} = q{} unless(exists $self->{_buffer_region_chr});
	unless(
	       $chr eq $self->{_buffer_region_chr} &&
	       $start_pos-$change_ref_offset >= $self->{_buffer_region_start} &&
	       $start_pos + $ref_seq_length + length($ref_right) < $self->{_buffer_region_end}
	       ){

		my $region_start = $start_pos-$change_ref_offset-10;
		my $region_end = $start_pos + $ref_seq_length + length($ref_right) + 5000;
		$self->{_buffer_region_chr} = $chr;
		$self->{_buffer_region_start} = $region_start;
		$self->{_buffer_region} = uc $fai->fetch($chr.':'.$region_start.'-'.$region_end); ##TODO why not work out the maximum amount of seq to pull back and just do it once - then just subsr?? Will save on IO...
		$self->{_buffer_region_end} = $region_start + length($self->{_buffer_region});
	}

	my $_buffer_region = $self->{_buffer_region};
	my $_buffer_region_start = $self->{_buffer_region_start};

	foreach my $read(@{$alignment}) {
		$read =~ s/ ([+-])/\t$1/ if($record_type eq 'D'); ## correction for a bug in the pindel output layout.... this is done here to allow use to pass a ref of the read into _parse_read
		_parse_read($record, $chr, $start_pos, \$read, $ref_seq_length, ($read_num++.$record_idx), $change_ref_offset, $change_ref_offset_end,\$_buffer_region,$_buffer_region_start, $self->{_noreads});
	}

	## This is not strictly read from the pindel input but is useful for woring out the number of repeats within the repeat-range.
	$record->repeats(_repeat_count($record,\$ref_left,\$ref_right));

	return 1;
}


=head3 _repeat_count

Takes the alignment string left and right components and attempts to resolve the number of repeats of the minimum change unit provided in the record.

@param1 = record    - a Sanger::CGP::Pindel::OutputGen::PindelRecord object. This must have the following fields
                      filled in: start, end, range_start, range_end, ref_seq, alt_seq, min_change.

@param2 = ref_left  - string_ref of the reference sequence before the change.

@param3 = ref_right - string_ref of the reference sequence after the change.

@return integer     - the number of times the min_change is seen within the range coordinates of the record.

=cut
sub _repeat_count {
	my ($record, $ref_left, $ref_right) = @_;

	my $min_change = $record->min_change || '';

	unless ($min_change) {
		print STDERR 'No minimum change found in record object.'."\n";
		print STDERR Dumper($record);
		croak 'No minimum change found in record object.'; ## jwh...
	}

	my $ref = $record->ref_seq || '';
	my $full_ref_string = $$ref_left . $ref . $$ref_right;
	my $relative_var_start = length($$ref_left);
	my $relative_var_end = $relative_var_start + length($ref);

	my $range_start_diff = $record->start - $record->range_start;
	my $range_end_diff   = $record->range_end - $record->end;

	my $relative_range_start = $relative_var_start - $range_start_diff + 1 ;
	my $relative_range_end = $relative_var_end + $range_end_diff - 1 ;

	my $repeat_range = substr($full_ref_string,$relative_range_start, ($relative_range_end - $relative_range_start));

	my $pre_rep_len = length $repeat_range;

	$repeat_range =~ s/($min_change)+//i;

	return int (($pre_rep_len-(length $repeat_range)) / length $min_change);
}

=head3 most_prev_change

Loops through an array of alignment strings sub-stringing the variant out using the
offsets provided. The most prevalent variant is returned. This was used when a bug
existed in the Pindel output where the event was selected from the first read in
the alignment list. This has since been corrected in newer versions.

@param1 = alignments - an array ref of alignment strings

@param2 = l_length - the left hand offset of where the variant begins

@param3 = c_length - the right hand offset of where the variant ends

@return (max_change_string, change_fraction) - the most prevalent change and its fraction of total reads at this loci

=cut
sub most_prev_change {
	my ($alignments, $l_length, $c_length) = @_;
	my %changes;
	my $total = 0;
	my $max_change = '';
	my $max_change_count = 0;

	foreach my $align(@{$alignments}) {
		my $tmp_change = substr($align, $l_length, $c_length);

		if(++$changes{$tmp_change} >= $max_change_count){
			$max_change = $tmp_change;
			$max_change_count = $changes{$tmp_change};
		}
		$total++;
	}

	return ($max_change, ($max_change_count / $total));
}

=head3 _shrink_change

Takes a string and attempts to reduce it to its repetitive component if it has one i.e.

ATCATCATC shrinks to ATC
ATCATCATCG will not shrink

Repeats upto $max_repeat_unit_size in length are looked for.

@param1 = change               - the change sequence to shrink.

@param2 = max_repeat_unit_size - the maximum repeat unit length.

@returns string - the change string reduced to its repetative component

=cut
sub _shrink_change {
	my ($change, $max_repeat_unit_size) = @_;
	my $new_change;
	my $change_sub_length = 0;
	my ($test_change, $chg_tmp);
	if(length $change > 1) {
		my $last_size = length $change;
		for my $sub_len(1..$max_repeat_unit_size) {
			$test_change = substr($change,0,$sub_len);
			$chg_tmp = $change;
			$chg_tmp =~ s/^($test_change)+//i;
			unless (length $chg_tmp) {
				$new_change = $test_change;
				last;
			}
		}
	}
	## return the shrunken change otherwise the original input change
	return defined $new_change ? $new_change : $change;
}

=head3 _parse_read

Takes a read string in the pidnel output and adds it to the record object as a sam record.
We will need these reads to do two things. 1) create pindel-bam files 2) count unique reads
between bwa and pindel.

@param1 record           - a Sanger::CGP::Pindel::OutputGen::PindelRecord object. This must have the following fields
                           filled in: type, idx, chro, start, ref_seq, alt_seq.

@param2 chro             - The chromosome name.

@param3 start_pos        - The genomic start position of the event.


@param4 read             - string_ref of a pindel read string.
                           <paddingSEQvarSEQpadding	strand	genomic-start-of-anchor-read	mapq-qual-of-anchor-read	sample-name	@read-name/read1|2>
                           (  TCCCCTaACAGTC   	+	1336638	37	COLO-829	@EAS188_62:3:72:801:1263/2 )

@param5 ref_seq_length   - length of the reference sequence associated with the event

@param6 read_idx         - the read identifier. This should be unique. This is used to form
                           part of the unique name of the read.

@param7 change_ref_start - this is the indicated start of the variant on the pindel record reference string and
                           is used to grab variant sequence from the read string.

@param8 change_ref_end   - this is the indicated end of the variant on the pindel record reference string and
                           is used to grab variant sequence from the read string.
=cut
sub _parse_read {
	my ($record, $chr, $start_pos, $read, $ref_seq_length, $read_idx, $change_ref_start, $change_ref_end, $_buffer_region, $_buffer_region_start, $no_read_data) = @_;
	$no_read_data ||= 0;

	my @bits = split /\t+/, ${$read};
	my $sample = $bits[-2];
	my $strand = $bits[-5];
	$strand =~ tr/\+\-/\-\+/; # need to invert as is strand of anchor

	my $read_data;
	if($no_read_data == 0) {
		## This is a custom read name component added to the read name when it is put into Pindel.
		## As pindel currently does not preserve the read group, if we want to identify read group
		## specific errors we need to track the read groups from the reads....
		my ($read_group) = $bits[-1] =~ /\/[12]_RG(.+)$/;
		$read_group = '' unless $read_group;

		my ($name, $rg_pair) = split /\//, $bits[-1];
		$name = substr($name,1) if substr($name,0,1) eq '@';

		# need this to force uniqness in reads that have multiple events
		# and make display in gbrowse work for overlapping reads
		$name .= '_r'.$read_idx;

		my $mapq = $bits[-3]; # mapq of anchor read, sensible to use this in the output

		# locate the left and right parts of the read string
		my $read_seq = $bits[0];
		my $read_left = substr($read_seq, 0, $change_ref_start);
		my $read_right = substr($read_seq, $change_ref_end);
		my $event = substr($read_seq, $change_ref_start,$change_ref_end - $change_ref_start);

		## we do this so that we can efficiently strip the space characters from the read components...
		my $left_seq_length = $read_left =~ tr/ATCGN/ATCGN/;## the tr simply counts the number of atcgs in the string... v-efficient
		my $right_seq_length = $read_right =~ tr/ATCGN/ATCGN/;## the tr simply counts the number of atcgs in the string... v-efficient
		my $event_seq_length = $event =~ tr/ATCGN/ATCGN/;
		my $event_length = length $event;

		## create a read sequence without any spaces. This is MUCH MUCH faster than using s///.
		my $space_stripped_read_seq = substr($read_left, (length($read_left) - $left_seq_length), $change_ref_start);
		$space_stripped_read_seq .= $event if $event_seq_length;
		$space_stripped_read_seq .= substr($read_right, 0,$right_seq_length);

		$start_pos -= $left_seq_length ;

		# Some data have -ve position starts so need to be corrected. i.e. the position starts before the beginning of the reference.
		# This occurs with species like Devil that map to shattered contig sequences.
		if($start_pos < 1) {
			my $corr_size = (abs $start_pos)+1;
			$read_seq = substr($space_stripped_read_seq, $corr_size);
			$read_left = substr(substr($read_left, (length($read_left) - $left_seq_length),$change_ref_start), $corr_size);

			$left_seq_length = $read_left =~ tr/ATCGN/ATCGN/;
			$start_pos = 1;
		}

	## Keep track of all the reads associated with a variant.
	## These bits are for pindel_bam creation.
	## These bam files only contain the reads identified from within pindel as having a variant
	## These bam files are used in things like gbrowse/jbrowse for display

		my @cig_list = ($left_seq_length, 'M');
		push @cig_list, $ref_seq_length, 'D' if $ref_seq_length;
		push @cig_list, $event_seq_length, 'I' if $event_seq_length;
		push @cig_list, $right_seq_length, 'M';

		my $flag = $strand eq '-' ? 16 : 0; # previously 1+8 but as not really paired anymore
		my @tags = calmd($chr, $start_pos, \@cig_list, \$space_stripped_read_seq, $_buffer_region, $_buffer_region_start);

		pop @tags; # we dont need the last value.. at the moment...
		unshift @tags, "RG:Z:$read_group" if($read_group);
		$read_data = [$name,$flag,$chr,$start_pos,$mapq,join(q{},@cig_list),'*','0','0',$space_stripped_read_seq,'*',@tags];
	}

    # Create a basic sam line for the read and add it to the record.
    $record->add_read($sample,$strand,$read_data);

	return 1;
}

=head3 calmd

Calculate the MD and NM tags based on new Cigar string, seq and genomic location
Added to deal with problem in samtools calmd where split reads are handled poorly

@param1 chr             - chromosome name.

@param2 start           - genomic start position of the read. If the read is soft clipped the start position should take this into account.

@param3 cigar           - array_ref containing cigar components of the read.

@param4 seq_ref         - a string_ref of the read sequence.

@param5 reference_ref   - a string_ref of the reference sequence. This sequence must cover the entirety of the read and any deletions.

@param6 reference_start - genomic start position of the reference sequence.

@returns                - (MD string, NM string).

=cut
sub calmd {
	my ($chr, $start, $cigar, $seq_ref, $reference_ref, $reference_start) = @_;

	my $pairs_length = scalar @{$cigar};

	my @md_bits;
	my $nm = 0;
	my $read_seq_idx = 0;
	my $was_previous_type_match = 0;
	my ($type, $len, $new_start, $ref_seq, $ref_base, $read_base, $match_count, $final_end);

	for(my $i=0;$i<$pairs_length;$i+=2){
		$len  = $cigar->[$i];
		$type = $cigar->[$i+1];

		if($type eq 'I') {
			$nm += $len;
			$read_seq_idx+=$len;
			#$was_previous_type_match = 0; ## I's do not show up in MD tag...
			next;
		}

		if($type eq 'S') {
			# does not affect NM
			$read_seq_idx+=$len;
			$was_previous_type_match = 0;
			next;
		}
		# start does not move when an insert or soft clip
		$new_start = $start + $len;

		if($type eq 'N') {## what the heck is this??? jwh this will never be created as part of the _parse_read method...
			$start = $new_start; # need to do this here as well
			$read_seq_idx+=$len;
			next;
		}
		$final_end = $new_start-1;

		#$ref_seq = $fai->fetch($chr.':'.$start.'-'.$final_end);
		## this should grab the event ref seq from the huge ref_string we are holding for this region.....
		## we then sub string each component from this ref_string... here we assume that the ref string will always be large enough as this should have been worked out elsewhere....
		## This may mean that we have to bring the function into an object method....
		$ref_seq = substr(${$reference_ref}, ($start-$reference_start), ($final_end - $start+1)) || '';

		$start = $new_start;
		if($type eq 'M') {
			$match_count = 0;
			foreach my $ref_seq_pos (0..length($ref_seq)-1) {
				$ref_base = substr($ref_seq,$ref_seq_pos,1);
				$read_base = substr(${$seq_ref},$read_seq_idx++,1);
				#warn "$ref_base, $read_base";
				if($read_base ne $ref_base){
					## If we get to this point we have found a single base pair difference between the reference and the M part of the read.
					## Mark it and carry on with the rest of the match length....
					$nm++;
					## stash the previous count...
					if($was_previous_type_match){
						$md_bits[-1] += $match_count; ## ok the previous type was a match so we need to merge them (skipping Ns?...);
					                                ## THIS BIT ALONE MEANS WE CANNOT SIMPLY CONCAT A STRING...
					}else{
						push @md_bits, $match_count;
					}

					$match_count = 0;
					push @md_bits, $ref_base, $match_count; ## This is a sub... start the next match region
					$was_previous_type_match = 1; ## The last thing in the md_bits array will now be a number so may need to append
				}else{
					$match_count++;
				}
			}

			## finally add the remaining match count to the md_bits array..
			if($was_previous_type_match){
				$md_bits[-1] += $match_count; ## ok the previous type was a match so we need to merge them (skipping Ns?...);
			}else{
				push @md_bits, $match_count;# if($match_count > 0);
				$was_previous_type_match = 1;
			}
			next;
		}

		if($type eq 'D') {
		  push @md_bits, '^'.$ref_seq;
			$nm += $len;
			$was_previous_type_match = 0;
			next;
		}

		die "Cigar contains unhandled type ($type): $cigar";
	}

	my $md = join q{}, @md_bits;
	return ('MD:Z:'.$md, 'NM:i:'.$nm, $final_end);
}


sub calmd_orig {
	my ($chr, $start, $cigar, $seq_ref, $fai) = @_;
	my @c_lengths = split /[^[:digit:]]/xms, $cigar;
	my @c_types = split /[[:digit:]]+/xms, $cigar;
	shift @c_types; # as always has empty first element

	my @md_bits;
	my $nm = 0;
	my ($type, $len, @ref, $new_start, $ref_seq, $ref_base, $read_base, $match_count, $final_end);
	my @read_seq = split //xms, ${$seq_ref};
	for my $e(0..((scalar @c_lengths)-1)) {
	$type = $c_types[$e];
	$len = $c_lengths[$e];

	if($type eq 'I') {
		$nm += $len;
		splice(@read_seq,0,$len); # remove seq from read
		next;
	}

	if($type eq 'S') {
		# does not affect NM
		splice(@read_seq,0,$len); # remove seq from read
		next;
	}
	# start does not move when an insert or soft clip
	$new_start = $start + $len;

	if($type eq 'N') {
		$start = $new_start; # need to do this here as well
		next;
	}
	$final_end = $new_start-1;
	$ref_seq = $fai->fetch($chr.':'.$start.'-'.$final_end); ##TODO why not work out the maximum amount of seq to pull back and just do it once - then just subsr?? Will save on IO...

	$start = $new_start;
	if($type eq 'M') {
		@ref = split //xms, $ref_seq;
		$match_count = 0;
		while($ref_base = shift @ref) {
			$read_base = shift @read_seq;
			if($read_base eq $ref_base) {
				$match_count++;
				next;
			}
			$nm++;
			push @md_bits, $match_count if($match_count > 0);
			push @md_bits, $ref_base;
			$match_count = 0;
		}
		push @md_bits, $match_count if($match_count > 0);
		next;
	}
	if($type eq 'D') {
		push @md_bits, '^'.$ref_seq;
		$nm += $len;
		next;
	}

	die "Cigar contains unhandled type ($type): $cigar";

	}
	my @final_md;
	foreach my $m_bit(@md_bits) {
		if(@final_md == 0) {
			push @final_md, 0 if($m_bit !~ m/^[0-9]+$/);
			push @final_md, $m_bit;
			next;
		}
		if($m_bit =~ m/^[0-9]+$/ && $final_md[-1] =~ m/^[0-9]+$/) {
			$final_md[-1] += $m_bit;
			next;
		}
		if($m_bit !~ m/^[0-9]+$/ && $final_md[-1] !~ m/^[0-9]+$/) {
			push @final_md, 0, $m_bit;
			next;
		}
		push @final_md, $m_bit;
	}
	push @final_md, 0 if($final_md[-1] !~ m/^[0-9]+$/);
	my $md = join q{}, @final_md;
	return ('MD:Z:'.$md, 'NM:i:'.$nm, $final_end);
}
