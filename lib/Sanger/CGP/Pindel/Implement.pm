package Sanger::CGP::Pindel::Implement;

use strict;
use warnings FATAL => 'all';
use autodie qw(:all);
use Const::Fast qw(const);
use File::Spec;
use File::Which qw(which);
use File::Path qw(make_path);
use File::Temp qw(tempfile);
use Capture::Tiny;
use List::Util qw(first);
use FindBin qw($Bin);

use Sanger::CGP::Pindel;
our $VERSION = Sanger::CGP::Pindel->VERSION;

use PCAP::Threaded;
use PCAP::Bam;

use Sanger::CGP::Pindel::OutputGen::BamUtil;

const my $PINDEL_GEN_COMM => ' -b %s -o %s -t %s';
const my $SAMTOOLS_FAIDX => ' faidx %s %s > %s';
const my $FILTER_PIN_COMM => ' %s %s %s %s';
const my $PINDEL_COMM => ' %s %s %s %s %s %s';
const my $PIN_2_VCF => q{ -mt %s -wt %s -r %s -d %s -o %s -so %s -mtp %s -wtp %s -pp '%s'};

my $l_bin = $Bin;

sub prepare {
  my $options = shift;
  $options->{'tumour_name'} = (PCAP::Bam::sample_name($options->{'tumour'}))[0];
  $options->{'normal_name'} = (PCAP::Bam::sample_name($options->{'normal'}))[0];
  return 1;
}

sub input {
  my ($index, $options) = @_;
  return 1 if(exists $options->{'index'} && $index != $options->{'index'});

  my $tmp = $options->{'tmp'};
  return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), $index);

  my @inputs = ($options->{'tumour'}, $options->{'normal'});
  my $iter = 1;
  for my $input(@inputs) {
    next if($iter++ != $index); # skip to the relevant input in the list

    ## build command for this index
    #

    my $max_threads = int ($options->{'threads'} / scalar @inputs);

    my $sample = sanitised_sample_from_bam($input);
    my $gen_out = File::Spec->catdir($tmp, $sample);
    make_path($gen_out);

    my $command = "$^X ";
    $command .= _which('pindel_input_gen.pl');
    $command .= sprintf $PINDEL_GEN_COMM, $input, $gen_out, $max_threads;

    PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, $index);

    #
    ## The rest is auto-magical
    PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), $index);
  }
  return 1;
}

sub split {
  my ($index, $options) = @_;
  return 1 if(exists $options->{'index'} && $index != $options->{'index'});

  my $tmp = $options->{'tmp'};
  return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), $index);

  my @seqs = sort keys %{$options->{'seqs'}};
  my $iter = 1;
  for my $seq(@seqs) {
    next if($iter++ != $index); # skip to the relevant seq in the list

    ## build command for this index
    #

    my $gen_out = File::Spec->catdir($tmp, 'refs');
    make_path($gen_out);

    my $command = _which('samtools');
    $command .= sprintf $SAMTOOLS_FAIDX,  $options->{'reference'},
                                          $seq,
                                          File::Spec->catfile($gen_out, "$seq.fa");

    PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, $index);

    #
    ## The rest is auto-magical
    PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), $index);
  }
  return 1;
}

sub filter {
  my ($index, $options) = @_;
  return 1 if(exists $options->{'index'} && $index != $options->{'index'});

  my $tmp = $options->{'tmp'};
  return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), $index);

  my @seqs = sort keys %{$options->{'seqs'}};
  my $iter = 1;
  for my $seq(@seqs) {
    next if($iter++ != $index); # skip to the relevant seq in the list

    ## build command for this index
    #

    my $gen_out = File::Spec->catdir($tmp, 'filter');
    make_path($gen_out);

    my $refs = File::Spec->catdir($tmp, 'refs');

    my $command = _which('filter_pindel_reads');
    $command .= sprintf $FILTER_PIN_COMM, File::Spec->catfile($refs, "$seq.fa"),
                                          $seq,
                                          File::Spec->catfile($gen_out, $seq),
                                          (join q{ }, @{$options->{'seqs'}->{$seq}});

    PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, $index);

    #
    ## The rest is auto-magical

    PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), $index);
  }
  return 1;
}

sub pindel {
  my ($index, $options) = @_;
  return 1 if(exists $options->{'index'} && $index != $options->{'index'});

  my $tmp = $options->{'tmp'};
  return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), $index);

  my @seqs = sort keys %{$options->{'seqs'}};
  my $iter = 1;
  for my $seq(@seqs) {
    next if($iter++ != $index); # skip to the relevant seq in the list

    ## build command for this index
    #

    my $gen_out = File::Spec->catdir($tmp, 'pout');
    make_path($gen_out);

    my $filtered = File::Spec->catdir($tmp, 'filter');
    my $refs = File::Spec->catdir($tmp, 'refs');
    my ($bd_fh, $bd_file) = tempfile('pindel_db_XXXX', UNLINK => 1);
    close $bd_fh;

    my $command = _which('pindel');
    $command .= sprintf $PINDEL_COMM, File::Spec->catfile($refs, "$seq.fa"),
                                      File::Spec->catfile($filtered, $seq),
                                      $gen_out,
                                      $seq,
                                      $bd_file,
                                      5;

    PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, $index);

    # a little cleanup
    for my $ext((qw(BP INV LI TD))) {
      unlink File::Spec->catfile($gen_out, (join '_', $seq, $seq, $ext));
    }
    unlink $bd_file;

    #
    ## The rest is auto-magical

    PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), $index);
  }
  return 1;
}

sub pindel_to_vcf {
  my $options = shift;

  my $tmp = $options->{'tmp'};
  return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 0);

  ## build command for this index
  #

  my $pg = Sanger::CGP::Pindel::OutputGen::BamUtil::pg_from_caller('pindel', 'cgpPindel indel detection', $VERSION, $options->{'cmd'});
  $pg =~ s/\t/\\t/g;

  my $command = $^X.' '._which('pindel_2_combined_vcf.pl');
  $command .= sprintf $PIN_2_VCF, $options->{'tumour'},
                                  $options->{'normal'},
                                  $options->{'reference'},
                                  File::Spec->catdir($tmp, 'pout'), # pindel result dir
                                  File::Spec->catfile($tmp, 'pindel.vcf'),
                                  File::Spec->catfile($tmp, 'pindel'),
                                  $options->{'seqtype'},
                                  $options->{'seqtype'},
                                  $pg;
  # optional items
  $command .= ' -s' if(defined $options->{'skipgerm'});
  $command .= ' -as '.$options->{'assembly'} if(defined $options->{'assembly'});
  $command .= ' -sp '.$options->{'species'}   if(defined $options->{'species'});

  PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, 0);

  # any cleanup


  #
  ## The rest is auto-magical

  PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 0);
}

sub sanitised_sample_from_bam {
  my $sample = (PCAP::Bam::sample_name(shift))[0];
  $sample =~ s/[^\/a-z0-9_-]/_/ig; # sanitise sample name
  return $sample;
}

sub determine_jobs {
  my $options = shift;
  my $tmp = $options->{'tmp'};
  my @exclude;
  @exclude = split /,/, $options->{'exclude'} if(exists $options->{'exclude'});
  my %seqs;
  for my $in_bam($options->{'tumour'}, $options->{'normal'}) {
    my $samp_path = File::Spec->catdir($tmp, sanitised_sample_from_bam($in_bam));
    my @files = file_list($samp_path, qr/\.txt$/);
    for my $file(@files) {
      my ($seq) = $file =~ m/(.+)\.txt$/;
      next if(first { $seq eq $_ } @exclude);
      push @{$seqs{$seq}}, File::Spec->catfile($samp_path, $file);
    }
  }
  $options->{'seqs'} = \%seqs;
  my @seqs = keys %seqs;
  return scalar @seqs;
}

sub file_list {
  my ($dir, $regex) = @_;
  my @files;
  opendir(my $dh, $dir) || die;
  while(readdir $dh) {
    push @files, $_ if($_ =~ $regex);
  }
  closedir $dh;
  return @files;
}

sub _which {
  my $prog = shift;
  my $path = File::Spec->catfile($l_bin, $prog);
  $path = which($prog) unless(-e $path);
  return $path;
}

1;
