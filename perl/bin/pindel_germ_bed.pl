#!/usr/bin/perl

########## LICENCE ##########
# Copyright (c) 2014 Genome Research Ltd.
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


BEGIN {
  use Cwd qw(abs_path);
  use File::Basename;
  unshift (@INC,dirname(abs_path($0)).'/../lib');
};

use strict;
use warnings FATAL => 'all';
use autodie qw(:all);

use File::Path qw(remove_tree make_path);
use Getopt::Long;
use Pod::Usage qw(pod2usage);
use Const::Fast qw(const);
use IO::Zlib;

use PCAP::Cli;
use Sanger::CGP::Pindel;
use Vcf;

const my $DEFAULT_FLAG => 'F012';

my $options = &setup;

my $o_fh;
if(exists $options->{'output'} && defined $options->{'output'}) {
  open $o_fh, '>', $options->{'output'};
}
else {
  $o_fh = *STDOUT;
}

my $flag = $options->{'flag'};

my $vcf = Vcf->new( file=>$options->{'input'},
                    version=>'4.1');
$vcf->parse_header();

my $has_entry = 0;
while(my $x = $vcf->next_data_array){
  next unless(";$x->[6];" =~ m/;$flag;/); # skip things that don't have this flag value
  my ($ref, $start, $alt) = (@{$x})[0,1,4];
  my $end = $start + length $alt;
  # start is always the base before the change in VCF in/del so this magically ends up as half-open coord.
  print $o_fh (join "\t", $ref, $start, $end),"\n" or die "Failed to write: $!";
  $has_entry = 1;
}
#Fake line as first position of first contig unless we have had at least one bed entry
if($has_entry == 0){
  my $first_contig = (sort keys %{$vcf->get_header_line(key=>'contig')->[0]})[0];
  print $o_fh (join ("\t",$first_contig,0,1),"\n") or die "Failed to write: $!";
}


close $o_fh;

sub setup {
  my %opts;
  pod2usage(-msg  => "\nERROR: Option must be defined.\n", -verbose => 1,  -output => \*STDERR) if(scalar @ARGV == 0);
  $opts{'cmd'} = join " ", $0, @ARGV;
  GetOptions( 'h|help' => \$opts{'h'},
              'm|man' => \$opts{'m'},
              'i|input=s' => \$opts{'input'},
              'o|output=s' => \$opts{'output'},
              'v|version' => \$opts{'version'},
              'f|flags=s' => \$opts{'flag'},

  ) or pod2usage(2);

  $opts{'flag'} = $DEFAULT_FLAG unless(defined $opts{'flag'});

  pod2usage(-verbose => 1) if(defined $opts{'h'});
  pod2usage(-verbose => 2) if(defined $opts{'m'});

  if($opts{'version'}) {
    print 'Version: ',Sanger::CGP::Pindel->VERSION,"\n";
    exit 0;
  }

  PCAP::Cli::file_for_reading('input', $opts{'input'});

  return \%opts;
}

__END__

=head1 pindel_germ_bed.pl

Generates a bed file of likely germline indels for this sample, used as a filter input for caveman.

=head1 SYNOPSIS

pindel_germ_bed.pl [options]

  Required parameters:
    -input     -i   Flagged Pindel VCF file

  Optional
    -output    -o   Path to save BED file to [stdout]
    -flags     -f   Filter flag that denotes a germline event [F012]

  Other:
    -help      -h   Brief help message.
    -man       -m   Full documentation.
    -version   -v   Version

  Example:
    pindel_germ_bed.pl -i cgpPindel.flagged.vcf -o germline.bed -f F012

