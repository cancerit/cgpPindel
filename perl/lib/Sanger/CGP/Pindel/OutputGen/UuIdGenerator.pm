package Sanger::CGP::Pindel::OutputGen::UuIdGenerator;

use Sanger::CGP::Pindel;
our $VERSION = Sanger::CGP::Pindel->VERSION;

use strict;
use Data::UUID;

1;

sub new{
	my $proto = shift;
	my (%args) = @_;
	my $class = ref($proto) || $proto;

	my $self = {
		_gen => new Data::UUID,
	};
    bless $self, $class;
    return $self;
}

sub next{
	my $gen = shift->{_gen};
	return lc $gen->to_string($gen->create);
}

sub reset{}
