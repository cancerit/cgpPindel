package Sanger::CGP::Pindel::InputGen::PairToPindel;

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
use Carp qw(croak);

use Sanger::CGP::Pindel;

const my $I_SIZE => 500;
const my $PINDEL_REC => qq{\@%s/%s_RG%s\n%s\n%s\t%s\t%s\t%s\t%s\t%s};

sub new {
  my ($class, $sample_name, $rg_pis) = @_;
  my $self = { 'sample' => $sample_name };
  bless $self, $class;
  $self->setup_rg_inserts($rg_pis) if(defined $rg_pis);
  return $self;
}

sub setup_rg_inserts {
  my ($self, $rg_pis) = @_;
  # this is purely to get round limitations of threads
  my %rgs;
  my @rg_sets = split /,/, $rg_pis;
  for (@rg_sets) {
    my ($rg, $pi) = $_ =~ m/^(.*):([[:digit:]]+)$/;
    $rgs{$rg} = $pi;
  }
  $self->{'rgs'} = \%rgs;
  return 1;
}

sub pair_to_pindel {
  my ($self, $pair_in) = @_;
  $self->{'pair'} = $pair_in;
  my $r1 = $self->{'pair'}->{'r1'};
  my $r2 = $self->{'pair'}->{'r2'};
  my @pindel_records;

  if($r1->good_query) {
    if($r2->good_anchor) {
      push @pindel_records, $self->_anchored_to_pindel($r2, $r1);
    }
    elsif($r1->good_anchor) {
      push @pindel_records, $self->_self_anchored_to_pindel($r1);
    }
  }
  if($r2->good_query) {
    if($r1->good_anchor) {
      push @pindel_records, $self->_anchored_to_pindel($r1, $r2);
    }
    elsif($r2->good_anchor) {
      push @pindel_records, $self->_self_anchored_to_pindel($r2);
    }
  }
  return \@pindel_records;
}

sub _anchored_to_pindel {
  my ($self, $anchor, $query) = @_;
  my $rg = $query->get('rg');
  my $isize = $self->{'rgs'}->{$rg};
  die "Failed to get insert size for readgroup: '$rg'\n" unless($isize);
  my $anchor_pos;
  if($anchor->strand eq '+') {
    $anchor_pos = $anchor->get('pos');
  }
  else {
    $anchor_pos = $anchor->get('pos') + $anchor->mapped_seq;
  }
  return sprintf $PINDEL_REC, $query->get('qname')
                            , $query->read_end
                            , $rg
                            , $query->seq_for_pindel
                            , $anchor->strand
                            , $anchor->get('rname')
                            , $anchor_pos # pindel says 3'
                            , $anchor->get('mapq')
                            , $isize
                            , $self->{'sample'};
}

sub _self_anchored_to_pindel {
  my ($self, $query) = @_;
  my $rg = $query->get('rg');
  my $isize = $self->{'rgs'}->{$rg};
  die "Failed to get insert size for readgroup: '$rg'\n" unless($isize);
  my $anchor_strand = ($query->strand eq '+') ? '-' : '+';
  my $anchor_pos;
  if($anchor_strand eq '+') {
#    $anchor_pos = $query->get('pos') - ($isize*2) + $query->mapped_seq;
    $anchor_pos = $query->get('pos') - $isize;# - length $query->get('seq');
  }
  else {
#    $anchor_pos = $query->get('pos') + ($isize*2) + $query->mapped_seq;
    $anchor_pos = $query->get('pos') + length $query->mapped_seq;
  }
  return sprintf $PINDEL_REC, $query->get('qname')
                            , $query->read_end
                            , $rg
                            , $query->seq_for_pindel
                            , $anchor_strand
                            , $query->get('rname')
                            , $anchor_pos  # pindel says 3'
                            , $query->get('mapq')
                            , $isize
                            , $self->{'sample'};
}

1;

__DATA__

