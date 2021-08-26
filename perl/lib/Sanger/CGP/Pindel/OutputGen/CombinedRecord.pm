# Copyright (c) 2014-2021
#
# Author: CASM/Cancer IT <cgphelp@sanger.ac.uk>
#
# This file is part of cgpPindel.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
# 1. The usage of a range of years within a copyright statement contained within
# this distribution should be interpreted as being equivalent to a list of years
# including the first and last year specified and all consecutive years between
# them. For example, a copyright statement that reads ‘Copyright (c) 2005, 2007-
# 2009, 2011-2012’ should be interpreted as being identical to a statement that
# reads ‘Copyright (c) 2005, 2007, 2008, 2009, 2011, 2012’ and a copyright
# statement that reads ‘Copyright (c) 2005-2012’ should be interpreted as being
# identical to a statement that reads ‘Copyright (c) 2005, 2006, 2007, 2008,
# 2009, 2010, 2011, 2012’.
#
package Sanger::CGP::Pindel::OutputGen::CombinedRecord;

use Sanger::CGP::Pindel;

use strict;
use base qw{Sanger::CGP::Pindel::OutputGen::PindelRecord};

1;

sub new{
	my $proto = shift;
	my (%args) = @_;
	my $class = ref($proto) || $proto;
	my $self = $class->SUPER::new(@_);

	$self->{_valid} = $args{'-valid'};

	$self->{_p_mt_pos} = $args{'-p_mt_pos'};
	$self->{_p_mt_neg} = $args{'-p_mt_neg'};
	$self->{_p_wt_pos} = $args{'-p_wt_pos'};
	$self->{_p_wt_neg} = $args{'-p_wt_neg'};

	$self->{_b_mt_pos} = $args{'-b_mt_pos'};
	$self->{_b_mt_neg} = $args{'-b_mt_neg'};
	$self->{_b_wt_pos} = $args{'-b_wt_pos'};
	$self->{_b_wt_neg} = $args{'-b_wt_neg'};

	$self->{_d_mt_pos} = $args{'-d_mt_pos'};
	$self->{_d_mt_neg} = $args{'-d_mt_neg'};
	$self->{_d_wt_pos} = $args{'-d_wt_pos'};
	$self->{_d_wt_neg} = $args{'-d_wt_neg'};

	$self->{_rd_mt_pos} = $args{'-rd_mt_pos'};
	$self->{_rd_mt_neg} = $args{'-rd_mt_neg'};
	$self->{_rd_wt_pos} = $args{'-rd_wt_pos'};
	$self->{_rd_wt_neg} = $args{'-rd_wt_neg'};

	$self->{_uc_mt_pos} = $args{'-uc_mt_pos'};
	$self->{_uc_mt_neg} = $args{'-uc_mt_neg'};
	$self->{_uc_wt_pos} = $args{'-uc_wt_pos'};
	$self->{_uc_wt_neg} = $args{'-uc_wt_neg'};

	$self->{_fd_wt} = $args{'-fd_wt'};
	$self->{_fc_wt} = $args{'-fc_wt'};
	$self->{_fd_mt} = $args{'-fd_mt'};
	$self->{_fc_mt} = $args{'-fc_mt'};

	return $self;
}

sub valid{
	my($self,$value) = @_;
	$self->{_valid} = $value if defined $value;
	return $self->{_valid};
}

sub generic_setter {
	my ($self, $target, $value) = @_;
	return $self->$target($value);
}

sub p_mt_pos{
	my($self,$value) = @_;
	$self->{_p_mt_pos} = $value if defined $value;
	return $self->{_p_mt_pos};
}

sub p_wt_pos{
	my($self,$value) = @_;
	$self->{_p_wt_pos} = $value if defined $value;
	return $self->{_p_wt_pos};
}

sub p_mt_neg{
	my($self,$value) = @_;
	$self->{_p_mt_neg} = $value if defined $value;
	return $self->{_p_mt_neg};
}

sub p_wt_neg{
	my($self,$value) = @_;
	$self->{_p_wt_neg} = $value if defined $value;
	return $self->{_p_wt_neg};
}

sub b_mt_pos{
	my($self,$value) = @_;
	$self->{_b_mt_pos} = $value if defined $value;
	return $self->{_b_mt_pos};
}

sub b_wt_pos{
	my($self,$value) = @_;
	$self->{_b_wt_pos} = $value if defined $value;
	return $self->{_b_wt_pos};
}

sub b_mt_neg{
	my($self,$value) = @_;
	$self->{_b_mt_neg} = $value if defined $value;
	return $self->{_b_mt_neg};
}

sub b_wt_neg{
	my($self,$value) = @_;
	$self->{_b_wt_neg} = $value if defined $value;
	return $self->{_b_wt_neg};
}

sub d_mt_pos{
	my($self,$value) = @_;
	$self->{_d_mt_pos} = $value if defined $value;
	return $self->{_d_mt_pos};
}

sub d_wt_pos{
	my($self,$value) = @_;
	$self->{_d_wt_pos} = $value if defined $value;
	return $self->{_d_wt_pos};
}

sub d_mt_neg{
	my($self,$value) = @_;
	$self->{_d_mt_neg} = $value if defined $value;
	return $self->{_d_mt_neg};
}

sub d_wt_neg{
	my($self,$value) = @_;
	$self->{_d_wt_neg} = $value if defined $value;
	return $self->{_d_wt_neg};
}

sub rd_mt_pos{
	my($self,$value) = @_;
	$self->{_rd_mt_pos} = $value if defined $value;
	return $self->{_rd_mt_pos};
}

sub rd_wt_pos{
	my($self,$value) = @_;
	$self->{_rd_wt_pos} = $value if defined $value;
	return $self->{_rd_wt_pos};
}

sub rd_mt_neg{
	my($self,$value) = @_;
	$self->{_rd_mt_neg} = $value if defined $value;
	return $self->{_rd_mt_neg};
}

sub rd_wt_neg{
	my($self,$value) = @_;
	$self->{_rd_wt_neg} = $value if defined $value;
	return $self->{_rd_wt_neg};
}

sub uc_mt_pos{
	my($self,$value) = @_;
	$self->{_uc_mt_pos} = $value if defined $value;
	return $self->{_uc_mt_pos};
}

sub uc_wt_pos{
	my($self,$value) = @_;
	$self->{_uc_wt_pos} = $value if defined $value;
	return $self->{_uc_wt_pos};
}

sub uc_mt_neg{
	my($self,$value) = @_;
	$self->{_uc_mt_neg} = $value if defined $value;
	return $self->{_uc_mt_neg};
}

sub uc_wt_neg{
	my($self,$value) = @_;
	$self->{_uc_wt_neg} = $value if defined $value;
	return $self->{_uc_wt_neg};
}

sub fd_wt{
	my($self,$value) = @_;
	$self->{_fd_wt} = $value if defined $value;
	return $self->{_fd_wt};
}

sub fd_mt{
	my($self,$value) = @_;
	$self->{_fd_mt} = $value if defined $value;
	return $self->{_fd_mt};
}

sub fc_wt{
	my($self,$value) = @_;
	$self->{_fc_wt} = $value if defined $value;
	return $self->{_fc_wt};
}

sub fc_mt{
	my($self,$value) = @_;
	$self->{_fc_mt} = $value if defined $value;
	return $self->{_fc_mt};
}
