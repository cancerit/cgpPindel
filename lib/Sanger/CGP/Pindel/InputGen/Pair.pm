package Sanger::CGP::Pindel::InputGen::Pair;

use strict;
use English qw( -no_match_vars );
use autodie qw(:all);
use warnings FATAL => 'all';
use Const::Fast qw(const);

use Sanger::CGP::Pindel::InputGen::Read;

const my $MIN_MAPQ => 0;
const my $SINGLE_END_MIN_MAPQ => 6;
const my $SOFTCLIP_MAX => 10;
const my $SOFTCLIP_MIN_AVG_PBQ => 10;

sub new {
  my ($class, $r1, $r2) = @_;
  my $self = {};
  bless $self, $class;
  $self->{'r1'} = Sanger::CGP::Pindel::Read->new($r1, 1);
  $self->{'r2'} = Sanger::CGP::Pindel::Read->new($r2, 2);
  return $self;
}

sub exact {
  my $self = shift;
  return 1 if($self->{'r1'}->exact && $self->{'r2'}->exact);
  return 0;
}

sub unmapped_pair {
  my $self = shift;
  return 1 if($self->{'r1'}->unmapped && $self->{'r2'}->unmapped);
  return 0;
}

sub has_good_anchor {
  my $self = shift;
  my $r2_state = $self->{'r2'}->good_anchor; # to ensure both fully populates
  return 1 if($self->{'r1'}->good_anchor || $r2_state);
  return 0;
}

sub has_good_query {
  my $self = shift;
  my $r2_state = $self->{'r2'}->good_query; # to ensure both fully populates
  return 1 if($self->{'r1'}->good_query || $r2_state);
  return 0;
}

sub keep_pair {
  my $self = shift;
  return 0 if($self->unmapped_pair);
  return 0 if($self->exact);
  return 0 unless($self->has_good_anchor);
  return 0 unless($self->has_good_query);
  return 1;
}

1;

__END__

=head1 NAME

Sanger::CGP::Pindel::Pair - Describes a pair of reads and their suitability as Pindel candidates

=head2 Constructor

=over 4

=item Sanger::CGP::Pindel::Pair->new($read1, $read2)

Create a pair object from a pair of sam lines

  my $pair = Sanger::CGP::Pindel::Pair->new($read1, $read2);
  next unless($pair->keep_pair);

=back

=head2 Methods

=over 4

=item parse_read

Convert SAM record into Sanger::CGP::Pindel::Read object.

=item exact

If both reads are exact mappings return true.

=item unmapped_pair

If both reads are unmapped return true.

=item has_good_anchor

Return true if one (or both) of the reads meets criteria of a good anchor.

=item has_good_query

Return true is one (or both) of the reads meets criteria of a good query.

=item keep_pair

Return true if read in pair are suitable for conversion to a Pindel Input Record

=back
