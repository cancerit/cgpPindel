package Sanger::CGP::Pindel::PairToPindel;

use strict;
use English qw( -no_match_vars );
use autodie qw(:all);
use warnings FATAL => 'all';
use Const::Fast qw(const);
use Carp qw(croak);

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
    my ($rg, $pi) = split /:/, $_;
    $rgs{$rg} = $pi;
  }
  $self->{'rgs'} = \%rgs;
  return 1;
}

sub set_pair {
  my ($self, $pair) = @_;
  croak "No pair provided" unless(defined $pair);
  $self->{'pair'} = $pair;
}

sub pair_to_pindel {
  my $self = shift;
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
  my $isize = $self->{'rgs'}->{$query->get('rg')};
  my $anchor_pos;
  if($anchor->strand eq '+') {
    $anchor_pos = $anchor->get('pos');
  }
  else {
    $anchor_pos = $anchor->get('pos') + $anchor->mapped_seq;
  }
  return sprintf $PINDEL_REC, $query->get('qname')
                            , $query->read_end
                            , $query->get('rg') || 0
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
  my $isize = $self->{'rgs'}->{$query->get('rg')};
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
                            , $query->get('rg') || 0
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

