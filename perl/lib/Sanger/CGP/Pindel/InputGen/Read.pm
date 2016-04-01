package Sanger::CGP::Pindel::InputGen::Read;

########## LICENCE ##########
# Copyright (c) 2014-2016 Genome Research Ltd.
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


use strict;
use English qw( -no_match_vars );
use autodie qw(:all);
use warnings FATAL => 'all';
use Const::Fast qw(const);
use List::Util qw(sum sum0);
use List::MoreUtils qw(first_index last_index);

use Sanger::CGP::Pindel;

# remove before release
use Carp qw(croak);
use Data::Dumper;

const my $PROPER_MAPPED => 2;
const my $UNMAPPED      => 4;
const my $READ_REVERSED => 16; #0x0010 (set when reverse)
const my $FIRST_IN_PAIR => 64;
const my $SECOND_IN_PAIR => 128;

const my $MIN_MAPQ => 0;
const my $MIN_ANCHOR_MAPQ => 0;
const my $MIN_ANCHOR_MAPPED => 21;

const my $POOR_QUAL_PBQ => 19;
const my $MAX_POOR_PBQ_FRAC => 0.4;

const my $MAX_CLIP_BOTH_END => 10;
const my $MAX_SOFTCLIP_FRAC => 0.7;
const my $MAX_CIGAR_OPS_FOR_ANCHOR => 7*2; #cigar operations array both elements


sub new {
  my ($class, $sam, $end, $tabix) = @_;
#  my @elements = (split /\t/, ${$sam})[0,1,2,3,4,5,8,9,10];
  my ($qname, $flag, $rname, $pos, $mapq, $cigar, $seq, $qual) = (split /\t/, ${$sam})[0,1,2,3,4,5,9,10];
  # just clean this up as it is of no use
  $cigar =~ s/[[:digit:]]+H//g if(index($cigar, 'H') != -1);
  my ($rg) = ${$sam} =~ m/\tRG:Z:([^\t]+)/;
  my $self = {'qname' => $qname,
              'flag' => int $flag,
              'rname' => $rname,
              'pos' => int $pos,
              'mapq' => int $mapq,
              'cigar' => $cigar,
              'seq' => $seq,
              'qual' => $qual,
              'end' => int $end,
              'rg' => defined $rg ? $rg : '.',
              };
  $self->{'tabix'} = $tabix if(defined $tabix);
  bless $self, $class;
  $self;
}

sub get {
  my ($self, $key) = @_;
  if(exists $self->{$key}) {
    return $self->{$key};
  }
  return undef;
}

sub read_end {
  my $self = shift;
  return 1 if(($self->{'flag'} | $FIRST_IN_PAIR) == $self->{'flag'});
  return 2 if(($self->{'flag'} | $SECOND_IN_PAIR) == $self->{'flag'});
  die Dumper($self)."\n\nERROR: Reads are unpaired acording to flags";
}

sub unmapped {
  my $self = shift;
  (($self->{'flag'} | $UNMAPPED) == $self->{'flag'}) ? 1 : 0;
}

sub reversed {
  my $self = shift;
  (($self->{'flag'} | $READ_REVERSED) == $self->{'flag'}) ? 1 : 0;
}

sub strand {
  my $self = shift;
  $self->reversed ? '-' : '+';
}

sub exact {
  my $self = shift;
  return $self->{'_exact'} if(exists $self->{'_exact'});
  if(($self->{'flag'} | $PROPER_MAPPED) != $self->{'flag'}) {
    $self->{'_exact'} = 0;
  }
  elsif($self->{'cigar'} !~ m/^[[:digit:]]+M$/) {
    $self->{'_exact'} = 0;
  }
  else {
    $self->{'_exact'} = 1;
  }
  $self->{'_exact'};
}

sub _softclip_cal {
  my $self = shift;
  return 1 if(exists $self->{'_softclip_leading'}); # only need 1 to be present to know calc has occurred
  my @cigar_operations = @{$self->_cigar_operations};

  $self->{'_softclip_leading'} = $cigar_operations[0] if($cigar_operations[1] eq 'S');
  $self->{'_softclip_trailing'} = $cigar_operations[-2] if($cigar_operations[-1] eq 'S');

  $self->{'_softclip_leading'} ||= 0;
  $self->{'_softclip_trailing'} ||= 0;

  1;
}

sub good_anchor {
  my $self = shift;
  return $self->{'_good_anchor'} if(exists $self->{'_good_anchor'});
  $self->{'_good_anchor'} = $self->_good_anchor;
  $self->{'_good_anchor'};
}

sub _good_anchor {
  my $self = shift;
  return 0 if($self->unmapped);
  return 0 if($self->{'mapq'} <= $MIN_ANCHOR_MAPQ);
  return 0 if($self->mapped_seq <= $MIN_ANCHOR_MAPPED);
  return 0 if((scalar @{$self->_cigar_operations}) > $MAX_CIGAR_OPS_FOR_ANCHOR);
  return 0 if($self->frac_pbq_poor > $MAX_POOR_PBQ_FRAC);
  return 0 if(exists $self->{'tabix'} && $self->_tabix_hit);
  1;
}

sub _tabix_hit {
  my $self = shift;
  my $tabix = $self->{'tabix'};

  my $iter = $self->{'tabix'}->query(sprintf '%s:%d-%d', $self->{'rname'}, $self->{'pos'}-1, $self->{'pos'});
  return 0 unless(defined $iter);
  while(my $ret = $iter->next){
    return 1;
  }
  return 0;
}

sub frac_pbq_poor {
  my $self = shift;
  unless(exists $self->{'_frac_pbq_poor'}) {
    my $quals = $self->full_qual_array;
    my $bad = 0;
    for(@{$quals}) {
      $bad++ if($_ <= $POOR_QUAL_PBQ);
    }
    $self->{'_frac_pbq_poor'} = $bad / (scalar @{$quals});
  }
  $self->{'_frac_pbq_poor'};
}

sub full_qual_array {
  my $self = shift;
  $self->{'_full_qual_array'} = [map{ord($_)- 33 } split //, $self->{'qual'}] unless(exists $self->{'_full_qual_array'});
  $self->{'_full_qual_array'};
}

sub good_query {
  my $self = shift;
  $self->{'_good_query'} = $self->_good_query unless(exists $self->{'_good_query'});
  $self->{'_good_query'};
}

sub _good_query {
  my $self = shift;
  return 0 if(index($self->{'seq'},'NN') >= 0);
  return 0 if($self->frac_pbq_poor > $MAX_POOR_PBQ_FRAC);
  unless($self->unmapped) {
    return 0 if($self->exact);
    return 0 if($self->softclip_frac > $MAX_SOFTCLIP_FRAC);
  }

  1;
}

# only use when deciding on anchor usage
# for example high softclip frac is bad if single end mapping
# BUT if single end and > 0.4 read is mapped can fake anchor for this read (combined with hard limit of seq)
sub softclip_frac {
  my $self = shift;
  unless(exists $self->{'_softclip_frac'}) {
    $self->{'_softclip_frac'} = $self->softclip_total / $self->mappable_seq;
  }
  $self->{'_softclip_frac'};
}

sub softclip_total {
  my $self = shift;
  unless(exists $self->{'_softclip_total'}) {
    $self->{'_softclip_total'} = $self->softclip_leading + $self->softclip_trailing;
  }
  $self->{'_softclip_total'};
}

sub softclip_leading {
  my $self = shift;
  $self->_softclip_cal unless(exists $self->{'_softclip_leading'});
  $self->{'_softclip_leading'};
}

sub softclip_trailing {
  my $self = shift;
  $self->_softclip_cal unless(exists $self->{'_softclip_trailing'});
  $self->{'_softclip_trailing'};
}

sub mappable_seq {
  my $self = shift;
  $self->{'_mappable_seq'} = length $self->{'seq'} unless(exists $self->{'_mappable_seq'});
  $self->{'_mappable_seq'};
}

sub mapped_seq {
  my $self = shift;
  $self->{'_mapped_seq'} = (length $self->{'seq'}) - $self->softclip_total unless(exists $self->{'_mapped_seq'});
  $self->{'_mapped_seq'};
}

sub qual_trim_seq {
  my $self = shift;
  my @full_qual = @{$self->full_qual_array};
  my $max = (scalar @full_qual) - 1;
  my $pos = 0;
  for (0..($max-1)) {
    if($full_qual[$_] > $POOR_QUAL_PBQ && $full_qual[$_+1] > $POOR_QUAL_PBQ) {
      $pos = $_;
      last;
    }
  }
  my $first_good = $pos;

  my @rev_qual = reverse @full_qual;
  $pos = 0;
  for (0..($max-1)) {
    if($rev_qual[$_] > $POOR_QUAL_PBQ && $rev_qual[$_+1] > $POOR_QUAL_PBQ) {
      $pos = $_;
      last;
    }
  }
  my $last_good = scalar @full_qual - $pos;

  return $self->{'seq'} if($first_good == 0 && $last_good + 1 == scalar @full_qual);
  my $capture = ($last_good - $first_good)+1;
  my $new_seq = substr($self->{'seq'}, $first_good, $capture);
  $new_seq;
}

sub seq_for_pindel {
  my $self = shift;
  unless(exists $self->{'_seq_for_pindel'}) {
    my $new_seq = $self->qual_trim_seq;
    if($self->reversed) {
      $new_seq =~ tr/ACGT/TGCA/;
      $new_seq = reverse $new_seq;
    }
    $self->{'_seq_for_pindel'} = $new_seq;
  }
  $self->{'_seq_for_pindel'};
}

sub _cigar_operations {
  my $self = shift;
  # need to write test for '=' as cigar op
  $self->{'_cigar_operations'} = [$self->{'cigar'} =~ m/([[:digit:]]+)([[:upper:]=])/g] unless(exists $self->{'_cigar_operations'});
  $self->{'_cigar_operations'};
}

1;
