package Sanger::CGP::Pindel::InputGen::Pair;

########## LICENCE ##########
# Copyright (c) 2014-2018 Genome Research Ltd.
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


use strict;
use English qw( -no_match_vars );
use autodie qw(:all);
use warnings FATAL => 'all';
use Const::Fast qw(const);

use Sanger::CGP::Pindel;

use Sanger::CGP::Pindel::InputGen::Read;

const my $MIN_MAPQ => 0;
const my $SINGLE_END_MIN_MAPQ => 6;
const my $SOFTCLIP_MAX => 10;

sub new {
  my ($class, $r1, $r2, $tabix) = @_;
  my $self = {};
  bless $self, $class;
  $self->{'r1'} = Sanger::CGP::Pindel::InputGen::Read->new($r1, 1, $tabix);
  $self->{'r2'} = Sanger::CGP::Pindel::InputGen::Read->new($r2, 2, $tabix);
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

sub qcfailed_pair {
  my $self = shift;
  return 1 if($self->{'r1'}->qc_failed && $self->{'r2'}->qc_failed);
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
  return 0 if($self->qcfailed_pair);
  return 0 if($self->unmapped_pair);
  return 0 if($self->exact);
  return 0 unless($self->has_good_anchor);
  return 0 unless($self->has_good_query);
  return 1;
}

1;

__END__

=head1 Sanger::CGP::Pindel::InputGen::Pair

Describes a pair of reads and their suitability as Pindel candidates

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
