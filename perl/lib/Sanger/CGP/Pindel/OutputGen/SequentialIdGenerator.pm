package Sanger::CGP::Pindel::OutputGen::SequentialIdGenerator;

use Sanger::CGP::Pindel;

use strict;

1;

sub new{
	my $proto = shift;
	my (%args) = @_;
	my $class = ref($proto) || $proto;

	my $self = {
		_counter => $args{'-start'} || 1,
	};
    bless $self, $class;
    return $self;
}

sub next{
	return shift->{_counter}++;
}

sub reset{
	shift->{_counter} = shift;
}
