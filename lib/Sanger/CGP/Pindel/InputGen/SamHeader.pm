package Sanger::CGP::Pindel::InputGen::SamHeader;

use strict;
use English qw( -no_match_vars );
use autodie qw(:all);
use warnings FATAL => 'all';
use Const::Fast qw(const);
use Carp qw( croak );

use Sanger::CGP::Pindel;
our $VERSION = Sanger::CGP::Pindel->VERSION;

use PCAP::Bam::Bas;

const my $DEFAULT_PI => 500;

sub new {
  my ($class, $bam) = @_;
  my $self = {'header_lines' => [],};
  $self->{'bas'} = PCAP::Bam::Bas->new("$bam.bas") if(defined $bam && -e "$bam.bas");
  bless $self, $class;
  return $self;
}

sub add_header_line {
  my ($self, $line) = @_;
  croak "add_header_line requires value for 'line'" unless(defined $line);
  chomp $line;
  croak "add_header_line received blank input" unless(length $line > 0);
  croak "add_header_line appears to have received invalid input $line" if($line !~ m/^\@/);
  return 1 if($line !~ m/^\@RG/);
  push @{$self->{'header_lines'}}, $line;
  return 1;
}

sub rgs_to_offsets_and_sample_name {
  my ($self, $force) = @_;
  if(scalar @{$self->{'header_lines'}} == 0) {
    die "No header information has been provided in the BAM file" unless($force);
    return {'.' => $DEFAULT_PI};
  }

  my %rgs;
  my $sample_name;
  for my $rg( @{$self->{'header_lines'}} ) {
    my ($rg_id) = $rg =~ m/\tID:([^\t]+)/;
    die "Readgroup line has no ID:\n\n$rg\n" unless(defined $rg_id);
    my ($median_insert) = $rg =~ m/\tPI:([[:digit:]]+)/;
    if(!defined $median_insert) {
      $median_insert = $self->{'bas'}->get($rg_id, 'median_insert_size') if(exists $self->{'bas'});
      if(!defined $median_insert) {
        warn "No PI tag found for RG $rg_id in header of *.bam.bas file";
        $median_insert = $DEFAULT_PI;
      }
    }
    die "BAM file has multiple RG lines with ID of: $rg_id" if(exists $rgs{$rg_id});
    $rgs{$rg_id} = $median_insert;

    my ($sample_tmp) = $rg =~ m/\SM:([^\t]+)/;
    if(defined $sample_name) {
      die "Different sample names found in RG lines, $sample_name vs. $sample_tmp" if($sample_name ne $sample_tmp);
    }
    else {
      $sample_name = $sample_tmp
    }
  }
  my $rg_str = q{};
  for (keys %rgs) {
    $rg_str .= "$_:$rgs{$_},";
  }
  chop $rg_str; # drop last comma.
  return ($rg_str, $sample_name);
}



1;

__END__

=head1 Sanger::CGP::Pindel::InputGen::SamHeader

Collate relevant header information.

=head2 Constructor

=over 4

=item Sanger::CGP::Pindel::SamHeader->new()

Create a SAM header parser

  my $sam_head = Sanger::CGP::Pindel::SamHeader->new();
  while(<>) {
    if($_ =~ m/^\@/) {
      $sam_head->add_header_line($_);
      next;
    }
    else {
      # this is just an example
      # you need more control code
      $sam_head->finalise;
    }
  }

=back

=head2 Methods

=over 4

=item add_header_line

Add a new header line to the object, will handle new line removal.

  $sam_head->add_header_line($line);

=item rgs_to_offsets

Takes the lines added and convert to useful data structure.

Warn on each readgroup that doesn't have a median insert size (PI) but use a default of 500.

Check for clashing RG IDs in header (should be rare).

=back
