package Sanger::CGP::Pindel::Implement;

########## LICENCE ##########
# Copyright (c) 2014-2018 Genome Research Ltd.
#
# Author: CASM/Cancer IT <cgphelp@sanger.ac.uk>
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
use warnings FATAL => 'all';
use autodie qw(:all);
use Const::Fast qw(const);
use File::Spec;
use File::Which qw(which);
use File::Copy qw(copy);
use File::Path qw(make_path remove_tree);
use File::Temp qw(tempfile);
use Capture::Tiny;
use List::Util qw(first);
use FindBin qw($Bin);

use Sanger::CGP::Pindel;

use PCAP::Threaded;
use PCAP::Bam;

use Sanger::CGP::Pindel::OutputGen::BamUtil;

const my $PINDEL_GEN_COMM => q{ -b %s -o %s -t %s};
const my $SAMTOOLS_FAIDX => q{ faidx %s %s > %s};
const my $FIFO_FILTER_TO_PIN => q{gunzip -c %s | %s %s %s %s /dev/stdin | %s %s %s %s %s %s %s};
const my $PIN_2_VCF => q{ -mt %s -wt %s -r %s -o %s -so %s -mtp %s -wtp %s -pp '%s' -i %s};
const my $PIN_MERGE => q{ -o %s -i %s -r %s};
const my $FLAG => q{ -a %s -u %s -s %s -i %s -o %s -r %s};
const my $PIN_GERM => q{ -f %s -i %s -o %s};
const my $BASE_GERM_RULE => 'F012'; # prefixed with additional F if fragment filtering.

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

    my $max_threads = $options->{'threads'};
    unless (exists $options->{'index'}) {
      $max_threads = int ($max_threads / scalar @inputs);
    }
    $max_threads = 1 if($max_threads == 0);

    my $sample = sanitised_sample_from_bam($input);
    my $gen_out = File::Spec->catdir($tmp, $sample);
    make_path($gen_out) unless(-e $gen_out);

    my $command = "$^X ";
    $command .= _which('pindel_input_gen.pl');
    $command .= sprintf $PINDEL_GEN_COMM, $input, $gen_out, $max_threads;
    $command .= " -r $options->{reference}";
    $command .= " -e $options->{badloci}" if(exists $options->{'badloci'});

    PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, $index);

    #
    ## The rest is auto-magical
    PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), $index);
  }
  return 1;
}

sub pindel {
  my ($index_in, $options) = @_;
  my $tmp = $options->{'tmp'};

  return 1 if(exists $options->{'index'} && $index_in != $options->{'index'});

  my @seqs = sort keys %{$options->{'seqs'}};
  my @indicies = limited_indicies($options, $index_in, scalar @seqs);
	for my $index(@indicies) {
    next if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), $index);
    my $seq = $seqs[$index-1];

    ## build commands for this index
    #

    my @command_set;

    # was split
    my $refs = File::Spec->catdir($tmp, 'refs');
    make_path($refs) unless(-e $refs);

    my $refseq_file = File::Spec->catfile($refs, "$seq.fa");

    my $split_comm = _which('samtools');
    $split_comm .= sprintf $SAMTOOLS_FAIDX,  $options->{'reference'},
                                          $seq,
                                          $refseq_file;

    push @command_set, $split_comm;

    # was filter
    my $filter_out = File::Spec->catdir($tmp, 'filter');
    make_path($filter_out) unless(-e $filter_out);

    my $filtered_seq = File::Spec->catfile($filter_out, $seq);

    # pindel
    my $gen_out = File::Spec->catdir($tmp, 'pout');
    make_path($gen_out) unless(-e $gen_out);

    my ($bd_fh, $bd_file) = tempfile(File::Spec->catfile($tmp, 'pindel_db_XXXX'), UNLINK => 0);
    close $bd_fh;

    unlink $filtered_seq if(-e $filtered_seq);

    ## FIFO madness
    push @command_set, 'mkfifo '.$filtered_seq; # shell for this
    push @command_set, sprintf $FIFO_FILTER_TO_PIN, (join q{ }, @{$options->{'seqs'}->{$seq}}),
                                                _which('filter_pindel_reads'),
                                                $refseq_file,
                                                $seq,
                                                $filtered_seq,
                                                _which('pindel'),
                                                $refseq_file,
                                                $filtered_seq,
                                                $gen_out,
                                                $seq,
                                                $bd_file,
                                                5;

    PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), \@command_set, $index);

    # a little cleanup
    for my $ext((qw(BP INV LI TD))) {
      unlink File::Spec->catfile($gen_out, (join '_', $seq, $seq, $ext));
    }
    unlink $bd_file;
    unlink $refseq_file;
    unlink $filtered_seq;

    #
    ## The rest is auto-magical

    PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), $index);
  }
  return 1;
}

sub pindel_to_vcf {
  my ($index_in, $options) = @_;
  my $tmp = $options->{'tmp'};

  return 1 if(exists $options->{'index'} && $index_in != $options->{'index'});

  my @seqs = sort keys %{$options->{'seqs'}};
  my @indicies = limited_indicies($options, $index_in, scalar @seqs);
	for my $index(@indicies) {
    next if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), $index);
    my $seq = $seqs[$index-1];
    my $pout = File::Spec->catdir($tmp, 'pout');

    my @in_files;
    for my $type(qw(D SI)) {
      my $in_file = File::Spec->catfile($pout, $seq.'_'.$seq.'_'.$type);
      push @in_files, $in_file if(-e $in_file && -f $in_file);
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
  my $outstub = File::Spec->catfile($options->{'outdir'}, $options->{'tumour_name'}.'_vs_'.$options->{'normal_name'});
  my $command = "$^X ";
  $command .= _which('pindel_merge_vcf_bam.pl');
  $command .= sprintf $PIN_MERGE, $outstub, $vcf, $options->{'reference'};
  if($options->{'tumour'} =~ m/\.cram$/) {
    $command .= ' -c';
  }
  elsif(-e $options->{'tumour'}.'.csi') {
    $command .= ' -s';
  }

  PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, 0);

  PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 0);
}

sub flag {
  my $options = shift;

  my $tmp = $options->{'tmp'};
  return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 0);

# FlagVcf.pl
# -r ~kr2/GitHub/cgpPindel/perl/rules/genomicRules.lst
# -sr ~kr2/GitHub/cgpPindel/perl/rules/softRules.lst
# -a /lustre/scratch112/sanger/cgppipe/nst_pipe/test_ref/human/37/e58/vagrent/codingexon_regions.indel.bed.gz
# -u /lustre/scratch112/sanger/kr2/pan_cancer_test_sets/pindel_np_gen/huge_file.gff3.gz
# -s /lustre/scratch112/sanger/cgppipe/nst_pipe/test_ref/human/37/gsm_reference_repeat.gff.gz
# -i pindel_farm/PD13371a_vs_PD13371b.vcf.gz
# -o pindel_farm/PD13371a_vs_PD13371b.flag_new_np.github.vcf

  my $stub = File::Spec->catfile($options->{'outdir'}, $options->{'tumour_name'}.'_vs_'.$options->{'normal_name'});

  my $new_vcf = "$stub.flagged.vcf";
  my $command = "$^X ";
  $command .= _which('FlagVcf.pl');
  $command .= sprintf $FLAG,  $options->{'genes'},
                              $options->{'unmatched'},
                              $options->{'simrep'},
                              "$stub.vcf.gz", # input
                              $new_vcf, # output,
                              $options->{'filters'};
  $command .= ' -sr '.$options->{'softfil'} if(exists $options->{'softfil'} && defined $options->{'softfil'});
  $command .= ' -p '.$options->{'apid'} if(exists $options->{'apid'} && defined $options->{'softapidfil'});

  my $vcf_gz = $new_vcf.'.gz';
  my $bgzip = _which('bgzip');
  $bgzip .= sprintf ' -c %s > %s', $new_vcf, $vcf_gz;

  my $tabix = _which('tabix');
  $tabix .= sprintf ' -p vcf %s', $vcf_gz;

  my @commands = ($command, $bgzip, $tabix);

  my $germ_rule = find_germline_rule($options);
  if(defined $germ_rule) {
    my $germ_bed = "$stub.germline.bed";
    my $germ = "$^X ";
    $germ .= _which('pindel_germ_bed.pl');
    $germ .= sprintf $PIN_GERM, find_germline_rule($options), $vcf_gz, $germ_bed;
    push @commands, $germ;
  }

  PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), \@commands, 0);

  unlink $new_vcf;

  PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 0);
}

sub find_germline_rule {
  my $options = shift;
  open my $ffh, '<', $options->{'filters'} or die $!;
  my $filter;
  while(<$ffh>) {
    if($_ =~ m/(F?$BASE_GERM_RULE)/) {
      $filter = $1;
      last;
    }
  }
  close $ffh;
  return $filter;
}

sub prepare {
  my $options = shift;
  $options->{'tumour_name'} = (PCAP::Bam::sample_name($options->{'tumour'}))[0];
  $options->{'normal_name'} = (PCAP::Bam::sample_name($options->{'normal'}))[0];
  return 1;
}

sub valid_seqs {
  my $options = shift;
  my @good_seqs;

  my @all_seqs;
  open my $FAI_IN, '<', $options->{'reference'}.'.fai';
  while(<$FAI_IN>) { push @all_seqs, (split /\t/, $_)[0]; }
  close $FAI_IN;

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

sub sanitised_sample_from_bam {
  my $sample = (PCAP::Bam::sample_name(shift))[0];
  $sample =~ s/[^a-z0-9_-]/_/ig; # sanitise sample name
  return $sample;
}

sub determine_jobs {
  my $options = shift;
  my $tmp = $options->{'tmp'};
  my @valid_seqs = valid_seqs($options);
  my %seqs;
  for my $in_bam($options->{'tumour'}, $options->{'normal'}) {
    my $samp_path = File::Spec->catdir($tmp, sanitised_sample_from_bam($in_bam));
    my @files = file_list($samp_path, qr/\.txt(:?\.gz)$/);
    for my $file(@files) {
      my ($seq) = $file =~ m/(.+)\.txt(:?\.gz)$/;
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

sub limited_indicies {
	my ($options, $index_in, $count) = @_;
  my @indicies;
  if(exists $options->{'limit'}) {
    # main script checks index is not greater than limit or < 1
	  my $base = $index_in;
	  while($base <= $count) {
	    push @indicies, $base;
	    $base += $options->{'limit'};
	  }
	}
	else {
	  push @indicies, $index_in;
	}
	return @indicies;
}

sub fragmented_files {
  my ($base_dir, $files, $hprefix, $outfile) = @_;
  my $file_count = scalar @{$files};
  return $files if($file_count <= 100);
  warn "Extreme number of files, slow merging of $file_count files required\n";
  my $grep_non_header = qq{zgrep -v '^$hprefix' %s >> %s};
  my $header_count;

  my $merged_file = $base_dir.$outfile;
  my $is_first = 1;
  for my $file(@{$files}) {
    if($is_first == 1) {
      copy($base_dir.$file, $merged_file);
      $is_first = 0;
      next;
    }
    my $cmd = sprintf $grep_non_header, $base_dir.$file, $merged_file;
    warn "Running: $cmd\n";
    system([0,1],$cmd); # autodie in use, exit code 1 returned when no lines returned by grep
  }
  return [$merged_file];
}

sub _which {
  my $prog = shift;
  my $l_bin = $Bin;
  my $path = File::Spec->catfile($l_bin, $prog);
  $path = which($prog) unless(-e $path);
  return $path;
}

1;

__END__

=head1 NAME

Sanger::CGP::Pindel::Implement - Implementation functions for use with the L<PCAP::Threaded|https://github.com/ICGC-TCGA-PanCancer/PCAP-core/blob/dev/lib/PCAP/Threaded.pm> model.

=head1 METHODS

These are split into 2 sections.

L<Threaded Methods|"Threaded Methods"> are those that use the PCAP::Threaded model to run external processes.

L<Helper Methods|"Helper Methods"> are those that are either used internally or can be called from the main program.

=head2 Threaded Methods

These methods either use or are intended to be called via the L<PCAP::Threaded|https://github.com/ICGC-TCGA-PanCancer/PCAP-core/blob/dev/lib/PCAP/Threaded.pm> model.

=head3 input

Generates input for L<pindel|"pindel"> step.
To be executed within the L<PCAP::Threaded|https://github.com/ICGC-TCGA-PanCancer/PCAP-core/blob/dev/lib/PCAP/Threaded.pm> model.

Expected to be run as multiple concurrent (or single series processes), one per input BAM (i.e. 2).

  Arg[1]      : Hashref
  Example     : my $threads = PCAP::Threaded->new($options->{'threads'});
                ...
                $threads->add_function('input', \&Sanger::CGP::Pindel::Implement::input);
                ...
                $threads->run(2, 'input', $options);

=head3 pindel

Runs pindel over input files generated by L<input|"input">.
To be executed within the L<PCAP::Threaded|https://github.com/ICGC-TCGA-PanCancer/PCAP-core/blob/dev/lib/PCAP/Threaded.pm> model.

Expected to be run as multiple concurrent (or single series processes), one per reference sequence (bar any excluded items).

  Arg[1]      : Hashref
  Example     : my $threads = PCAP::Threaded->new($options->{'threads'});
                ...
                $threads->add_function('pindel', \&Sanger::CGP::Pindel::Implement::pindel);
                ...
                my $jobs = Sanger::CGP::Pindel::Implement::determine_jobs($options);
                ...
                $threads->run($jobs, 'pindel', $options);

=head3 pindel_to_vcf

Converts raw pindel text alignment output into VCF and SAM files.
To be executed within the L<PCAP::Threaded|https://github.com/ICGC-TCGA-PanCancer/PCAP-core/blob/dev/lib/PCAP/Threaded.pm> model.

Expected to be run as multiple concurrent (or single series processes), one per reference sequence (bar any excluded items).

  Arg[1]      : Hashref
  Example     : my $threads = PCAP::Threaded->new($options->{'threads'});
                ...
                $threads->add_function('pin2vcf', \&Sanger::CGP::Pindel::Implement::pindel_to_vcf);
                ...
                my $jobs = Sanger::CGP::Pindel::Implement::determine_jobs($options);
                ...
                $threads->run($jobs, 'pin2vcf', $options);


=head3 merge_and_bam

Builds command to merge VCF and BAM files from individual reference sequenced into single files using bin/pindel_merge_vcf_bam.html,
executed within the L<PCAP::Threaded|https://github.com/ICGC-TCGA-PanCancer/PCAP-core/blob/dev/lib/PCAP/Threaded.pm> model.

There should only be a single use of this method in a processing flow, normally the last step.

If this process has been marked as complete (existence of C<tmp/progress/Sanger::CGP::Pindel::Implement::merge_and_bam.0>) it will not execute.

On completion progress file is created (C<tmp/progress/Sanger::CGP::Pindel::Implement::merge_and_bam.0>)

  Arg[1]      : Hashref
  Example     : Sanger::CGP::Pindel::Implement::merge_and_bam(
                                                    {'tmp'         => $tmp_dir,
                                                     'tumour_name' => $name_from_tumour_bam,
                                                     'normal_name' => $name_from_normal_bam,
                                                     'outdir'      => $dir_for_final_output
                                                    });

=head2 Helper Methods

These methods are general purpose.

=head3 prepare

Extracts sample name from respective bam headers and adds to incoming hashref with key suffixed '_name'.

  Arg[1]      : Hashref
  Example     : prepare({'tumour'=> $tum_bam,
                        'normal'=> $norm_bam,});
  ReturnType  : 1 on success.
  Exceptions  : See PCAP::Bam::sample_name

L<PCAP::Bam::sample_name|https://github.com/ICGC-TCGA-PanCancer/PCAP-core/blob/dev/lib/PCAP/Bam.pm>

=head3 valid_seqs

Returns the list of reference sequences after removing any indicated in the exclude list.

  Arg[1]      : Hashref
  Example     : my @seqs = valid_seqs({'exclude' => $csv_seq_to_exclude,
                                       'reference' => $path_to_fasta_reference});

C<exclude> can include wildcard of C<%> which will match equivalently to C<.*> in a regex.

my $options = shift;

=head3 sanitised_sample_from_bam

Gets sample name from BAM header and replaces potentially problematic characters with '_'.

  Arg[1]      : Path to BAM file
  Example     : sanitised_sample_from_bam($bam_file_path);
  ReturnType  : Scalar (string)
  Exceptions  : See PCAP::Bam::sample_name

L<PCAP::Bam::sample_name|https://github.com/ICGC-TCGA-PanCancer/PCAP-core/blob/dev/lib/PCAP/Bam.pm>

=head3 determine_jobs

Calculate the number of jobs required for processing following pindel input generation taking into account
any sequences that should be ignored.  Specifically for paired tumour/normal pindel processing.

  Arg[1]      : Hashref
  Example     : my $jobs_required = determine_jobs({'tmp' => $tmp_dir,
                                                    'tumour' => $tumour_bam,
                                                    'normal' => $normal_bam,
                                                    'exclude' => $csv_seq_to_exclude,
                                                    'reference' => $path_to_fasta_reference});
  ReturnType  : Scalar (int)

=head3 file_list

Simple file collation without using expensive glob'ing.

  Arg[1]      : Path to existing directory
  Arg[2]      : regex to match files
  Example     : file_list($samp_path, qr/\.txt$/);
  ReturnType  : Array of file names (path not included)
  Exceptions  : Autodie based errors on failure to [open|close]dir
