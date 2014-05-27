#!/usr/bin/perl

BEGIN {
  use Cwd qw(abs_path);
  use File::Basename;
  unshift (@INC,dirname(abs_path($0)).'/../lib');
};

use strict;
use warnings FATAL => 'all';
use autodie qw(:all);
use Getopt::Long;
use Pod::Usage qw(pod2usage);
use Const::Fast qw(const);
use File::Which qw(which);
use Vcf;
use Sanger::CGP::Pindel;
use PCAP::Cli;


const my $PINDEL_RAW_MIN => 3;
const my $MIN_CALL_FRAC => 0.05;
const my $GFF_HEAD => "##gff-version 3\n##cgpPindel-version: %s\n";
const my $GFF_FORMAT => "%s\tpindel_np_from_vcf\tindel\t%s\t%s\t.\t.\tSAMPLE_COUNT=%s,%s\n";

use Data::Dumper;

my $options = &get_options;
my $calls = collate_calls($options);
my $gff3_file = output_gff($options, $calls);
bgzip_tabix($gff3_file);

sub collate_calls {
  my $options = shift;
  my $all_calls = {};
  for my $file (@{$options->{'files'}}) {
    my ($new_sample, $new_calls) = parse_vcf($options, $file);
    condense_calls($all_calls, $new_sample, $new_calls);
  }
  return $all_calls;
}

sub condense_calls {
  my ($all_calls, $new_sample, $new_calls) = @_;
  my $new = 0;
  for my $chr(sort keys %{$new_calls}) {
    for my $pos(sort {$a<=>$b} keys %{$new_calls->{$chr}}) {
      my ($final_calls, $final_depth) = (0,0);
      for my $counts(@{$new_calls->{$chr}->{$pos}}) {
        my ($uniq_calls, $uniq_depth) = @{$counts};
        $final_calls += $uniq_calls; # as pindel doesn't share reads between events (in this version)
        $final_depth = $uniq_depth if($uniq_depth > $final_depth); # as the largest event starting here will have the highest depth
      }
      $new++ unless(exists $all_calls->{$chr}->{$pos});
      push @{$all_calls->{$chr}->{$pos}}, "$new_sample=$final_calls|$final_depth";
    }
  }
  warn "New change starts: $new\t($new_sample)\n";
  return 1;
}

sub parse_vcf {
  my ($options, $file) = @_;
  my $calls = {};
  my $samp_id = $options->{'samp_id'};

  my $vcf = Vcf->new(file=>$file);
  $vcf->parse_header();
  my $req_sample = ($vcf->get_header_line(key=>'SAMPLE', ID=>$samp_id))[0][0]{'SampleName'};
  $vcf->set_samples(include=>[$samp_id]);
  my ($total, $kept, $new) = (0,0, 0);
  while(my $d = $vcf->next_data_hash) {
    $total++;
    my $gtypes = $d->{'gtypes'}->{$samp_id};

    next if(($gtypes->{'PP'} + $gtypes->{'NP'}) < $PINDEL_RAW_MIN);
    my $uniq_calls = $gtypes->{'PU'} + $gtypes->{'NU'};
    my $uniq_depth = $gtypes->{'PR'} + $gtypes->{'NR'};
    $uniq_depth ||= 1;
    next if($uniq_calls / $uniq_depth < $MIN_CALL_FRAC);
    my $gff_corrected_pos = $d->{'POS'} + 1;
    push @{$calls->{$d->{'CHROM'}}->{$gff_corrected_pos}}, [$uniq_calls,$uniq_depth];
    $kept++;
  }
  warn "Retained: $kept/$total\n";
  return ($req_sample, $calls);
}

sub output_gff {
  my ($options, $calls_in) = @_;
  my %calls = %{$calls_in}; # mainly for readability

  my $raw_gff3 = "$options->{output}.gff3";
  open my $gff3, '>', $raw_gff3;
  print $gff3 sprintf $GFF_HEAD, Sanger::CGP::Pindel->VERSION or die "Failed to write to $raw_gff3";
  print $gff3 "##$options->{cmd}\n" or die "Failed to write to $raw_gff3";
  for my $chr(sort keys %calls) {
    for my $pos(sort {$a<=>$b} keys %{$calls{$chr}}) {
      my @sample_sets = sort @{$calls{$chr}{$pos}};
      print $gff3 sprintf $GFF_FORMAT, $chr, $pos, $pos, scalar @sample_sets, (join q{,}, @sample_sets) or die "Failed to write to $raw_gff3";
    }
  }
  close $gff3;
  return $raw_gff3;
}

sub bgzip_tabix {
  my $gff3_file = shift;

  my $gff3_gz = "$gff3_file.gz";
  my $tbi = "$gff3_gz.tbi";
  unlink $gff3_gz if(-e $gff3_gz);
  unlink $tbi if(-e $tbi);

  # autodie makes this nice and clean
  my $bgzip = which('bgzip');
  system($bgzip, $gff3_file);

  my $tabix = which('tabix');
  system($tabix, '-p', 'gff', $gff3_gz);
  return 1;
}

sub get_options {
  my %opts;
  $opts{'cmd'} = join " ", $0, @ARGV;
  GetOptions( 'h|help'      => \$opts{'h'},
              'm|man'       => \$opts{'m'},
              'o|output=s'  => \$opts{'output'},
              's|samp_id=s' => \$opts{'samp_id'},
  ) or pod2usage(2);

  pod2usage(-verbose => 1) if(defined $opts{'h'});
  pod2usage(-verbose => 2) if(defined $opts{'m'});

  pod2usage(-message => q{ERROR: 'output' must be defined.}, -verbose => 1) unless(defined $opts{'output'});
  pod2usage(-message => q{ERROR: 'samp_id' must be defined.}, -verbose => 1) unless(defined $opts{'samp_id'});

  $opts{'files'} = \@ARGV;

  for my $file(@{$opts{'files'}}) {
    PCAP::Cli::file_for_reading($file, $file);
  }

  return \%opts;
}

##FORMAT=<ID=PR,Number=1,Type=Integer,Description="Total mapped reads on the positive strand">
##FORMAT=<ID=NR,Number=1,Type=Integer,Description="Total mapped reads on the negative strand">
##FORMAT=<ID=PU,Number=1,Type=Integer,Description="Unique calls on the positive strand">
##FORMAT=<ID=NU,Number=1,Type=Integer,Description="Unique calls on the negative strand">
##SAMPLE_COUNT=2,PD4975b=55|148,PD6722b=69|178

__END__

=head1 pindel_np_from_vcf.pl

Generate a pinel normal panel gff3 from VCF inputs

=head1 SYNOPSIS

pindel_np_from_vcf.pl [options] *.vcf[.gz]

  Required options:

    -output   -o  File stub to write final result to
                   - .gff3.gz will be appended.
    -samp_id  -s  Which sample ID should be used
                   - output from pindel.pl tags all data
                     as TUMOUR or NORMAL

  Other:
    -help     -h  Brief help message.
    -man      -m  Full documentation.
