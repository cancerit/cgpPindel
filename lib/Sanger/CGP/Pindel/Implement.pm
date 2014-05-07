package Sanger::CGP::Pindel::Implement;

use strict;
use warnings FATAL => 'all';
use autodie qw(:all);
use Const::Fast qw(const);
use File::Spec;
use File::Which qw(which);
use File::Path qw(make_path remove_tree);
use File::Temp qw(tempfile);
use Capture::Tiny;
use List::Util qw(first);
use FindBin qw($Bin);
use Capture::Tiny qw(capture_stdout);

use Sanger::CGP::Pindel;
our $VERSION = Sanger::CGP::Pindel->VERSION;

use PCAP::Threaded;
use PCAP::Bam;

use Sanger::CGP::Pindel::OutputGen::BamUtil;

const my $PINDEL_GEN_COMM => ' -b %s -o %s -t %s';
const my $SAMTOOLS_FAIDX => ' faidx %s %s > %s';
const my $FILTER_PIN_COMM => ' %s %s %s %s';
const my $PINDEL_COMM => ' %s %s %s %s %s %s';
const my $PIN_2_VCF => q{ -mt %s -wt %s -r %s -o %s -so %s -mtp %s -wtp %s -pp '%s' -i %s};
const my $PIN_MERGE => q{ -o %s %s %s %s};

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
    make_path($gen_out) unless(-e $gen_out);

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

sub valid_seqs {
  my $options = shift;
  my @good_seqs;
  my $fai_seqs = capture_stdout { system('cut', '-f', 1, $options->{'reference'}.'.fai' ); };
  my @all_seqs = split /\n/, $fai_seqs;
  if(exists $options->{'exclude'}) {
    my @exclude = split /,/, $options->{'exclude'};
    my @exclude_patt;
    for my $ex(@exclude) {
      $ex =~ s/%/.+/;
      push @exclude_patt, $ex;
    }

    for my $sq(@all_seqs) {
      push @good_seqs, $sq unless(first { $sq =~ m/^$_$/ } @exclude_patt);
    }
  }
  else {
    push @good_seqs, @all_seqs;
  }
  return @good_seqs;
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

    ## build commands for this index
    #

    my @command_set;

    # was split
    my $refs = File::Spec->catdir($tmp, 'refs');
    make_path($refs) unless(-e $refs);

    my $split_comm = _which('samtools');
    $split_comm .= sprintf $SAMTOOLS_FAIDX,  $options->{'reference'},
                                          $seq,
                                          File::Spec->catfile($refs, "$seq.fa");

    push @command_set, $split_comm;

    # was filter
    my $filter_out = File::Spec->catdir($tmp, 'filter');
    make_path($filter_out) unless(-e $filter_out);

    my $filter_comm = _which('filter_pindel_reads');
    $filter_comm .= sprintf $FILTER_PIN_COMM, File::Spec->catfile($refs, "$seq.fa"),
                                          $seq,
                                          File::Spec->catfile($filter_out, $seq),
                                          (join q{ }, @{$options->{'seqs'}->{$seq}});

    push @command_set, $filter_comm;

    # pindel
    my $gen_out = File::Spec->catdir($tmp, 'pout');
    make_path($gen_out) unless(-e $gen_out);

    my $filtered = File::Spec->catdir($tmp, 'filter');
    my ($bd_fh, $bd_file) = tempfile('pindel_db_XXXX', UNLINK => 1);
    close $bd_fh;

    my $pindel_comm = _which('pindel');
    $pindel_comm .= sprintf $PINDEL_COMM, File::Spec->catfile($refs, "$seq.fa"),
                                      File::Spec->catfile($filtered, $seq),
                                      $gen_out,
                                      $seq,
                                      $bd_file,
                                      5;

    push @command_set, $pindel_comm;

    PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), \@command_set, $index);

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
  my ($index, $options) = @_;
  return 1 if(exists $options->{'index'} && $index != $options->{'index'});
  my $tmp = $options->{'tmp'};
  return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), $index);

  my @seqs = sort keys %{$options->{'seqs'}};
  my $iter = 1;
  for my $seq(@seqs) {
    next if($iter++ != $index); # skip to the relevant seq in the list

    my $pout = File::Spec->catdir($tmp, 'pout');

    my @in_files;
    for my $type(qw(D SI)) {
      my $in_file = File::Spec->catfile($pout, $seq.'_'.$seq.'_'.$type);
      push @in_files, $in_file if(-s $in_file);
    }
    if(scalar @in_files > 0) {
      my $vcf = File::Spec->catdir($tmp, 'vcf');
      make_path($vcf) unless(-e $vcf);

      my $pg = Sanger::CGP::Pindel::OutputGen::BamUtil::pg_from_caller('pindel', 'cgpPindel indel detection', $VERSION, $options->{'cmd'});
      $pg =~ s/\t/\\t/g;

      my $command = $^X.' '._which('pindel_2_combined_vcf.pl');
      $command .= sprintf $PIN_2_VCF, $options->{'tumour'},
                                      $options->{'normal'},
                                      $options->{'reference'},
                                      File::Spec->catfile($vcf, $seq.'_pindel.vcf'),
                                      File::Spec->catfile($vcf, $seq.'_pindel'),
                                      $options->{'seqtype'},
                                      $options->{'seqtype'},
                                      $pg,
                                      (join ' -i ', @in_files);

      # optional items
      $command .= ' -s' if(defined $options->{'skipgerm'});
      $command .= ' -as '.$options->{'assembly'} if(defined $options->{'assembly'});
      $command .= ' -sp '.$options->{'species'}   if(defined $options->{'species'});

      PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, $index);
    }

    #
    ## The rest is auto-magical

    PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), $index);
  }
  return 1;
}

sub merge_and_bam {
  my $options = shift;

  my $tmp = $options->{'tmp'};
  return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 0);

  my $vcf = File::Spec->catdir($tmp, 'vcf');
  my $search_vcf = File::Spec->catfile($vcf, '*_pindel.vcf');
  my $search_mt = File::Spec->catfile($vcf, '*_pindel_mt.sam');
  my $search_wt = File::Spec->catfile($vcf, '*_pindel_wt.sam');
  my $outstub = File::Spec->catfile($options->{'outdir'}, $options->{'tumour_name'}.'_vs_'.$options->{'normal_name'});
  my $command = "$^X ";
  $command .= _which('pindel_merge_vcf_bam.pl');
  $command .= sprintf $PIN_MERGE, $outstub, $search_vcf, $search_mt, $search_wt;

  PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, 0);

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
  my @valid_seqs = valid_seqs($options);
  my %seqs;
  for my $in_bam($options->{'tumour'}, $options->{'normal'}) {
    my $samp_path = File::Spec->catdir($tmp, sanitised_sample_from_bam($in_bam));
    my @files = file_list($samp_path, qr/\.txt$/);
    for my $file(@files) {
      my ($seq) = $file =~ m/(.+)\.txt$/;
      if(first { $seq eq $_ } @valid_seqs) {
        push @{$seqs{$seq}}, File::Spec->catfile($samp_path, $file);
      }
    }
  }
  $options->{'seqs'} = \%seqs;
  my @seqs = keys %seqs;
  return scalar @seqs;
}

sub file_list {
  my ($dir, $regex) = @_;
  my @files;
  opendir(my $dh, $dir);
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
