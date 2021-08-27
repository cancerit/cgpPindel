#!/usr/bin/env perl
# Copyright (c) 2014-2021 Genome Research Ltd
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

use strict;
use warnings FATAL => 'all';
use autodie qw(:all);
use Cwd qw(abs_path);
use Pod::Usage qw(pod2usage);
use FindBin qw($Bin);
use lib "$Bin/../lib";
use Getopt::Long;
use Capture::Tiny qw(capture);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use Data::UUID;
use Set::IntervalTree;

=head
# Build the new header
zgrep -B 1000000 -m 1 '^#CHROM' old_vcfs/PD37237b.pindel.vcf.gz | head -n -1 > new_header
zgrep -hm 1 '^##SAMPLE' old_vcfs/PD37237b*.pindel.vcf.gz >> new_header
zgrep -m 1 '^#CHROM' old_vcfs/PD37237b.pindel.vcf.gz >> new_header
=cut

my $options = setup();
my ($vcf_by_sample, $sample_head) = vcf_by_samples($options->{vcfs});
my ($records, $sample_order) = collate_data($vcf_by_sample);

# write stuff
header($options->{output}, $options->{vcfs}, $sample_head);
records($options->{output}, $records, $sample_order, $options->{min_vaf}, $options->{np}, $options->{control}, $options->{min_blat});

sub records {
  my ($output, $records, $sample_order, $min_vaf, $np_tree, $control, $min_blat) = @_;
  my %ds = %{$records};
  my $uuid_gen = Data::UUID->new;
  my @samples = @{$sample_order};
  for my $chr(sort keys %ds) {
    my $chr_tree;
    if(defined $np_tree && exists $np_tree->{$chr}) {
      $chr_tree = $np_tree->{$chr};
    }
    else {
      $chr_tree = Set::IntervalTree->new();
    }
    for my $pos(sort {$a <=> $b} keys %{$ds{$chr}}) {
      for my $seq_key(sort keys %{$ds{$chr}{$pos}}) {
        next if($control && exists $ds{$chr}{$pos}{$seq_key}{$control});

        my ($info, $format) = @{$ds{$chr}{$pos}{$seq_key}{_DATA_}};

        if(defined $np_tree ) {
          $_ = q{;}.$info;
          my ($rs) = $_ =~ m/;RS=(\d+)/;
          my ($re) = $_ =~ m/;RE=(\d+)/;
          next if(@{$chr_tree->fetch($rs, $re)} > 0);
        }

        my ($ref, $alt) = split ':', $seq_key;
        my $row = join "\t", $chr, $pos, $uuid_gen->to_string($uuid_gen->create), $ref, $alt, q{.}, q{}, $info, $format;
        my $samples_with_min_vaf = 0;
        my $samples_with_min_blat = 0;
        for my $s(@samples) {
          if(exists $ds{$chr}{$pos}{$seq_key}{$s}) {
            #GT:S1:S2:PP:NP:WTP:WTN:WTM:MTP:MTN:MTM:VAF
            #./.:12:205.534:3:2:4:3:.:3:2:0.003:0.417
            my ($mtp, $mtn, $vaf) = (split /:/, $ds{$chr}{$pos}{$seq_key}{$s})[8,9,11];
            $row .= "\t".$ds{$chr}{$pos}{$seq_key}{$s};
            # need to review if these are all q{.} when one is
            $vaf = 0 if($vaf eq q{.});
            $mtp = 0 if($mtp eq q{.});
            $mtn = 0 if($mtn eq q{.});
            $samples_with_min_vaf++ if($vaf >= $min_vaf);
            $samples_with_min_blat++ if($mtp >= $min_blat && $mtn >= $min_blat);
          }
          else {
            $row .= "\t.";
          }
        }
        if($samples_with_min_vaf > 0 && $samples_with_min_blat > 0) {
          print $output $row."\n";
        }
      }
    }
  }
}

sub collate_data {
  my ($vcf_by_sample) = @_;
  my %ds;
  my @samples = sort keys %{$vcf_by_sample};
  for my $s(@samples) {
    my $fh = IO::Uncompress::Gunzip->new($vcf_by_sample->{$s}, MultiStream => 1, AutoClose=> 1) or die "gunzip failed: $GunzipError\n";;
    while (<$fh>) {
      next if(m/^#/);
      chomp;
      my ($chr, $pos, undef, $ref, $alt, undef, undef, $info, $format, $data) = split /\t/, $_;
      my $seq_key = sprintf '%s:%s', $ref, $alt;
      if(! exists $ds{$chr}{$pos}{$seq_key}) {
        $ds{$chr}{$pos}{$seq_key}{_DATA_} = [$info, $format];
      }
      $ds{$chr}{$pos}{$seq_key}{$s} = $data;
    }
  }
  return (\%ds, \@samples);
}

sub header {
  my ($output, $vcfs, $sample_head) = @_;
  my $metadata_header = _command_output(sprintf q{zgrep -B 1000000 -m 1 '^#CHROM' %s | grep -v '^##SAMPLE=' | head -n -1}, $vcfs->[0]);
  my $col_header = _command_output(sprintf q{zgrep -m 1 '^#CHROM' %s}, $vcfs->[0]);
  # remove last col as has sample will be added back
  ${$col_header} =~ s/\t[^\t]+$//;
  my @s_keys = sort keys %{$sample_head};

  print $output ${$metadata_header}."\n";
  for (@s_keys) {
    print $output ${$sample_head->{$_}}."\n";
  }
  print $output join "\t", ${$col_header}, @s_keys;
  print $output "\n";

  return 1;
}

sub vcf_by_samples {
  my ($vcfs) = @_;
  my %vcfs;
  my %samps;
  for my $v(@{$vcfs}) {
    my $samp_head = _command_output(sprintf q{zgrep -hm 1 '^##SAMPLE' %s}, $v);
    my ($sample) = ${$samp_head} =~ m/ID=([^,]+)/;
    $vcfs{$sample} = $v;
    $samps{$sample} = $samp_head;
  }
  return (\%vcfs, \%samps);
}

sub _command_output {
  my $command = shift;
  my ($c_out, $c_err, $c_exit) = capture { system($command); };
  if($c_exit) {
    warn "An error occurred while executing $command\n";
    warn "\tERROR$c_err\n";
    exit $c_exit;
  }
  chomp $c_out;
  return \$c_out;
}

sub np_lookup {
  my ($gff3, $min_samp) = @_;
  my $interval_count = 0;
  printf STDERR "Loading normal panel...\n";
  my %tree;
  my $z = IO::Uncompress::Gunzip->new($gff3, MultiStream => 1, AutoClose=> 1) or die "gunzip failed: $GunzipError\n";
  my $value = 1;
  while(my $line = <$z>) {
    next if ($line =~ m/^#/);
    chomp $line;
    # simple hash look up as only check start coord.
    my ($chr, $start, $info) = (split /\t/, $line)[0,3,7];
    my ($samp_count) = $info =~ m/^SAMPLE_COUNT=(\d+)/;
    next if($samp_count < $min_samp);

    $tree{$chr} = Set::IntervalTree->new() unless(exists $tree{$chr});
    $tree{$chr}->insert(\$value, $start, $start+1); # as half-open (i.e. last value is 1 past end)
    $interval_count++;
  }
  printf STDERR "\tdone, $interval_count intervals loaded.\n";
  return \%tree;
}

sub vcf_list {
  my $list_file = shift;
  my @vcfs;
  open my $LIST, '<', $list_file;
  map { chomp $_; push @vcfs, $_; } <$LIST>;
  close $LIST;
  return \@vcfs;
}

sub setup{
  my %opts = (
    'cmd' => join(" ", $0, @ARGV),
    'mnps' => 1,
    'min_vaf' => 0,
    'min_blat' => 4,
  );
  GetOptions( 'h|help' => \$opts{h},
              'm|man' => \$opts{m},
              'v|version' => \$opts{v},
              'o|output=s' => \$opts{output},
              'n|np=s' => \$opts{np},
              's|mnps=i' => \$opts{mnps},
              'k|min:f' => \$opts{min_vaf},
              'b|blat:i' => \$opts{min_blat},
              'd|debug' => \$opts{debug},
              'c|control=s' => \$opts{control},
              'l|list:s' => \$opts{list},
  );

  if(defined $opts{'v'}) {
    printf "Version: %s\n", Sanger::CGP::Pindel::Implement->VERSION;
    exit;
  }

  pod2usage(-verbose => 1) if(defined $opts{h});
  pod2usage(-verbose => 2) if(defined $opts{m});

  if(defined $opts{min_vaf}) {
    if($opts{min_vaf} < 0 || $opts{min_vaf} > 1) {
      print STDERR "ERROR: -vaf is a fraction and must be between 0 and 1.\n";
      pod2usage(-verbose => 1);
    }
    $opts{min_vaf} = (sprintf '%.3f', $opts{min_vaf}) + 0; # force number to be stored
  }

  my @vcfs;
  if($opts{list}) {
    @vcfs = @{vcf_list($opts{list})};
  }
  else {
    @vcfs = @ARGV;
  }
  $opts{vcfs} = \@vcfs;

  if(@vcfs < 2) {
    print STDERR "ERROR: More than 1 VCF input is required\n";
    pod2usage(-verbose => 1);
  }

  if(defined $opts{np}) {
    $opts{np} = np_lookup($opts{np}, $opts{mnps});
  }
  else {
    delete $opts{np};
  }

  $opts{align} = $opts{output}.'.sam' unless(defined $opts{align});

  open my $ofh, '>', $opts{output};
  $opts{output} = $ofh;

  return \%opts;
}

__END__

=head1 NAME

pindelCohortMerge.pl - Takes outputs from pindelCohort.pl and merges into a single file, filtering events.

=head1 SYNOPSIS

pindelCohortMerge.pl [options] A.vcf.gz B.vcf.gz...
or
pindelCohortMerge.pl [options] -l vcf.list

  Required parameters:
    -output    -o   File path for VCF output (not compressed)

  Optional parameters:
    -list      -l   Text list of vcf files, 1 per line
    -min       -k   Keep events VAF >= VALUE (3dp) for 1 or more samples
                     - default is to retain events even if VAF == 0/. for all samples
    -np        -n   Normal panel gff3 file - omit if no filtering required.
    -mnps      -s   Minimum normal panel samples required to exclude [default: 1]
    -control   -c   Exclude events where this sample has calls.
    -blat      -b   Exclude events where no sample has MTP >= N && MTN >= N [default: 4]

  Other:
    -help      -h   Brief help message.
    -man       -m   Full documentation.
    -version   -v   Prints the version number.

=head1 DESCRIPTION

B<pindelCohortMerge.pl> generate a vcf with collated samples.

For each common loci/REF/ALT merge the samples.

Not vcf-merge (from vcftools) as that fails to retain GT ordering, even between rows.

=cut
