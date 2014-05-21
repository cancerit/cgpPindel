package Sanger::CGP::Pindel::OutputGen::PindelRecord;

use Sanger::CGP::Pindel;

use strict;

1;

sub new{
	my $proto = shift;
	my (%args) = @_;
	my $class = ref($proto) || $proto;

	my $self = {
		-version => $args{'-version'},
		_id => $args{'-id'},
		_idx => $args{'-idx'},
		_reads => $args{'-reads'},
		_chro => $args{'-chro'},
		_start => $args{'-start'},
		_end => $args{'-end'},
		_range_start => $args{'-range_start'},
		_range_end => $args{'-range_end'},
		_length => $args{'-length'},
		_min_change => $args{'-min_change'},
		_lub => $args{'-lub'},
		_ref_seq => $args{'-ref_seq'},
		_alt_seq => $args{'-alt_seq'},
		_sum_ms => $args{'-sum_ms'},
		_s1 => $args{'-s1'},
		_s2 => $args{'-s2'},
		_type => $args{'-type'},
		_repeats => $args{'-repeats'},
		_num_samples => $args{'-num_samples'},
		_sample_contrib => $args{'-sample_contrib'},
	};
    bless $self, $class;
    return $self;
}

sub version{
	my($self,$value) = @_;
	$self->{_version} = $value if defined $value;
	return $self->{_version};
}

sub id{
	my($self,$value) = @_;
	$self->{_id} = $value if defined $value;
	return $self->{_id};
}

sub idx{
	my($self,$value) = @_;
	$self->{_idx} = $value if defined $value;
	return $self->{_idx};
}

=head reads
These are the reads collected form the pindel output that contribute to a call
A hash-ref in the form {sample(i.e. the tumour/normal)}->{strand(-/+)}->[[sam_record]].
e.g.
reads = {
	'PD12345' => {
		'-' => [[record1], [record2], [record3]....],
		'+' => [[record4], [record5], [record6]....]
	},
	'PD67890' => {
		'-' => [[record7], [record8], [record9]....],
		'+' => [[record10], [record11], [record12]....]
	}
}
=cut
sub reads{
	my($self,$value) = @_;
	$self->{_reads} = $value if defined $value;
	return $self->{_reads};
}

=head add_read
Adds a sam record for a contributing read.

@param1 $sample_name - the name of the sample associated with the read.

@param1 $strand - either '+' or '-' for forward or reverse strand.

@param1 $sam_record - a sam record array-ref.

=cut
sub add_read{
	my($self,$sample_name,$strand,$sam_record) = @_;
	push(@{$self->{_reads}->{$sample_name}->{$strand}}, $sam_record);
}

=head add_read
Returns an array-ref of sam records for a given sample and strand.

@param1 $sample_name - the name of the sample associated with the read.

@param1 $strand - either '+' or '-' for forward or reverse strand.

=cut
sub get_reads{
	my($self,$sample_name,$strand) = @_;
	return $self->{_reads}->{$sample_name}->{$strand};
}

=head samples
Returns an array of UNORDERED sample names asscociated with the record.

=cut
sub samples{
	my($self) = @_;
	return keys %{$self->{_reads}};
}

sub chro{
	my($self,$value) = @_;
	$self->{_chro} = $value if defined $value;
	return $self->{_chro};
}

sub start{
	my($self,$value) = @_;
	$self->{_start} = $value if defined $value;
	return $self->{_start};
}

sub end{
	my($self,$value) = @_;
	$self->{_end} = $value if defined $value;
	return $self->{_end};
}

sub range_start{
	my($self,$value) = @_;
	$self->{_range_start} = $value if defined $value;
	return $self->{_range_start};
}

sub range_end{
	my($self,$value) = @_;
	$self->{_range_end} = $value if defined $value;
	return $self->{_range_end};
}

sub length{
	my($self,$value) = @_;
	$self->{_length} = $value if defined $value;
	return $self->{_length};
}

sub min_change{
	my($self,$value) = @_;
	$self->{_min_change} = $value if defined $value;
	return $self->{_min_change};
}

sub lub{
	my($self,$value) = @_;
	$self->{_lub} = $value if defined $value;
	return $self->{_lub};
}

sub ref_seq{
	my($self,$value) = @_;
	$self->{_ref_seq} = $value if defined $value;
	return $self->{_ref_seq};
}

sub alt_seq{
	my($self,$value) = @_;
	$self->{_alt_seq} = $value if defined $value;
	return $self->{_alt_seq};
}

sub sum_ms{
	my($self,$value) = @_;
	$self->{_sum_ms} = $value if defined $value;
	return $self->{_sum_ms};
}

sub num_samples {
	my($self,$value) = @_;
	$self->{_num_samples} = $value if defined $value;
	return $self->{_num_samples};
}

sub sample_contrib {
	my($self,$value) = @_;
	$self->{_sample_contrib} = $value if defined $value;
	return $self->{_sample_contrib};
}

sub s1{
	my($self,$value) = @_;
	$self->{_s1} = $value if defined $value;
	return $self->{_s1};
}

sub s2{
	my($self,$value) = @_;
	$self->{_s2} = $value if defined $value;
	return $self->{_s2};
}

sub type{
	my($self,$value) = @_;
	$self->{_type} = $value if defined $value;
	return $self->{_type};
}

sub repeats{
	my($self,$value) = @_;
	$self->{_repeats} = $value if defined $value;
	return $self->{_repeats};
}
