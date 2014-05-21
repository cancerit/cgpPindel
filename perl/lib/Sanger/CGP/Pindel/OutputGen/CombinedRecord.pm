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

	$self->{_call_wt_rg_count} = $args{'-call_wt_rg_count'};
	$self->{_total_wt_rg_count} = $args{'-total_wt_rg_count'};
	$self->{_call_mt_rg_count} = $args{'-call_mt_rg_count'};
	$self->{_total_mt_rg_count} = $args{'-total_mt_rg_count'};

    return $self;
}

sub valid{
	my($self,$value) = @_;
	$self->{_valid} = $value if defined $value;
	return $self->{_valid};
}

sub read_groups{
	my($self,$value) = @_;
	$self->{_read_groups} = $value if defined $value;
	return $self->{_read_groups};
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

sub call_wt_rg_count{
	my($self,$value) = @_;
	$self->{_call_wt_rg_count} = $value if defined $value;
	return $self->{_call_wt_rg_count};
}

sub total_wt_rg_count{
	my($self,$value) = @_;
	$self->{_total_wt_rg_count} = $value if defined $value;
	return $self->{_total_wt_rg_count};
}

sub call_mt_rg_count{
	my($self,$value) = @_;
	$self->{_call_mt_rg_count} = $value if defined $value;
	return $self->{_call_mt_rg_count};
}

sub total_mt_rg_count{
	my($self,$value) = @_;
	$self->{_total_mt_rg_count} = $value if defined $value;
	return $self->{_total_mt_rg_count};
}
