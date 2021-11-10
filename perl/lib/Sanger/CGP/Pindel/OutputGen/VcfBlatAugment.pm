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
package Sanger::CGP::Pindel::OutputGen::VcfBlatAugment;
use strict;
use warnings FATAL => 'all';
use autodie qw(:all);
use Capture::Tiny qw(capture);
use Const::Fast qw(const);
use File::Basename;
use File::Path qw(remove_tree);
use File::Spec::Functions;
use File::Temp qw(tempfile tempdir);
use IO::Compress::Gzip qw(:constants gzip $GzipError);
use List::Util qw(min max);

use Bio::DB::HTS;
use Bio::DB::HTS::Faidx;
use Vcf;

use PCAP::Bam::Bas;
use Sanger::CGP::Pindel;
use Sanger::CGP::PindelPostProcessing::VcfSoftFlagger;
use Sanger::CGP::Vcf::VcfProcessLog;
use Sanger::CGP::Vcf::VcfUtil;

const my $SD_MULT => 2;
const my $V_FMT => 8;
const my $V_GT_START => 9;
const my $READS_ONLY => q{bash -c 'set -o pipefail; samtools view -uF 3840 %s %s | samtools fasta - > %s'};
const my $BLAT_ONLY => q{bash -c 'blat -t=dna -q=dna -noTrimA -minIdentity=95 -noHead -out=psl %s %s %s && pslPretty -long -axt %s %s %s /dev/stdout'};
# returns +ve and then -ve results
const my $LOCI_FMT => '%s:%d-%d';
const my $PAD_EVENT => 3;
const my $LARGE_D => 100;
const my $MAPPED_RL_MULT => 0.9;

1;

sub new{
  my $proto = shift;
  my (%args) = @_;
  my $class = ref($proto) || $proto;

  my $self = {
    input => $args{input},
    ref => $args{ref},
    ofh => $args{ofh}, # vcf
    hts_files => $args{hts_files},
    fill_in => $args{fill_in},
    debug => $args{debug} || 0,
  };
  bless $self, $class;

  $self->_init($args{outpath});

  return $self;
}

sub _init {
  my ($self, $outpath) = @_;
  my $vcf = Vcf->new(file => $self->{input});
  $vcf->parse_header;
  $self->{vcf} = $vcf;
  $self->_sample_order; # needs vcf
  $self->_align_output($outpath);
  $self->_add_headers unless($self->{fill_in});
  $self->_hts; # has to go before _buffer_sizes
  $self->_buffer_sizes; # has to go before _validate_sample
  $self->_validate_samples;
  # load the fai
  $self->{fai} = Bio::DB::HTS::Faidx->new($self->{ref});
  return 1;
}

sub _close_sams {
  my $self = shift;
  for my $sample(@{$self->{vcf_sample_order}}) {
    close $self->{sfh}->{$sample};
  }
  return 1;
}

sub _fa_dict {
  my $self = shift;
  open my $D, '<', $self->{ref}.'.dict' or die "Failed to find $self->{ref}.dict, please generate with 'samtools dict'";
  chomp(my @dict = <$D>);
  close $D;
  $self->{fa_dict} = \@dict;
}

sub _align_output {
  my ($self, $outpath) = @_;
  $self->_fa_dict();
  for my $sample(@{$self->{vcf_sample_order}}) {
    my $sam_file = catfile($outpath, (sprintf '%s.sam.gz', $sample));
    $self->{samfile}->{$sample} = $sam_file;
    $self->{bamfile}->{$sample} = catfile($outpath, (sprintf '%s.bam', $sample));
    unlink $sam_file if(-e $sam_file);

    my $SAM = new IO::Compress::Gzip $sam_file, -Level => Z_BEST_SPEED or die "IO::Compress::Gzip failed: $GzipError\n";
    print $SAM join "\n", @{$self->{fa_dict}};
    print $SAM "\n";
    $self->{sfh}->{$sample} = $SAM;
  }
  return 1;
}

sub _sample_order {
  my $self = shift;
  my $i = 9; # genotype sample col from 9
  my %samp_pos;
  my @ordered_samples = $self->{vcf}->get_samples;
  for my $s(@ordered_samples) {
    $samp_pos{$s} = $i++;
  }
  $self->{vcf_sample_pos} = \%samp_pos;
  $self->{vcf_sample_order} = \@ordered_samples;
  $self->{vcf_sample_count} = scalar @ordered_samples;
  return 1;
}

sub _sample_from_hts {
  my ($self, $hts) = @_;
  my $hts_sample;
  foreach my $line (split(/\n/,$hts->header->text)) {
    next unless($line =~ m/^\@RG/);
    chomp $line;
    ($hts_sample) = $line =~ m/SM:([^\t]+)/;
    last if(defined $hts_sample);
  }
  die sprintf "ERROR: Failed to find a SM tag in a readgroup header of %s\n", $hts->hts_path unless(defined $hts_sample);
  return $hts_sample;
}

sub _validate_samples {
  my $self = shift;
  # check BAM/CRAM and VCF have same sample
  my @samples = $self->{vcf}->get_samples();
  for my $vcf_s(sort @samples) {
    next if(exists $self->{hts}->{$vcf_s});
    die sprintf "ERROR: Sample '%s' is not represented in the BAM/CRAM files provided.\n", $vcf_s;
  }
  return 1;
}

sub _buffer_sizes {
  my $self = shift;
  my $max_ins = 0;
  my $max_rl = 0;
  my $min_rl = 1_000_000;
  for my $hts_sample(keys %{$self->{hts}}) {
    my $b = PCAP::Bam::Bas->new($self->{hts}->{$hts_sample}->hts_path.'.bas');
    my $sample;
    for my $rg($b->read_groups) {
      #my $m_sd = int ($b->get($rg, 'mean_insert_size') + ($b->get($rg, 'insert_size_sd') * $SD_MULT));
      my $m_sd = int ($b->get($rg, 'mean_insert_size') * 5);
      $max_ins = $m_sd if($m_sd > $max_ins);
      my $tmp_max = max ($b->get($rg, 'read_length_r1'), $b->get($rg, 'read_length_r2'));
      my $tmp_min = min ($b->get($rg, 'read_length_r1'), $b->get($rg, 'read_length_r2'));
      $min_rl = $tmp_min if($tmp_min < $min_rl);
      $max_rl = $tmp_max if($tmp_max > $max_rl);
      my $s = $b->get($rg, 'sample');
      if($sample) {
        die "ERROR: Multiple samples found in %s.bas\n", $self->{hts}->{$hts_sample}->hts_path if($sample ne $s);
        if($sample ne $hts_sample) {
          die "ERROR: Sample in bas file (%s) doesn't match bam/cram file (%s),
              %s vs %s.bas\n", $sample, $hts_sample,
              $self->{hts}->{$hts_sample}->hts_path, $self->{hts}->{$hts_sample}->hts_path;
        }
      }
      else {
        $sample = $s;
      }
    }
  }
  $self->{min_rl} = $min_rl;
  $self->{max_insert} = $max_ins;
  #$self->{target_pad} = $max_rl;
  $self->{target_pad} = $max_ins; # as expanded search space we need to expan the match space
  return 1;
}

sub _hts {
  my $self = shift;
  for my $hts(@{$self->{hts_files}}) {
    my $tmp = Bio::DB::HTS->new(-bam => $hts, -fasta => $self->{ref});
    my $sample = $self->_sample_from_hts($tmp);
    if(exists $self->{hts}->{$sample}) {
      die sprintf "ERROR: More than one BAM/CRAM file for sample %s, %s vs %s\n", $sample, $self->{hts}->{$sample}->hts_path, $hts;
    }
    $self->{hts}->{$sample} = $tmp;
  }
  return 1;
}

sub output_header {
  my $self = shift;
  my $fh = $self->{ofh};
  print $fh Sanger::CGP::PindelPostProcessing::VcfSoftFlagger::reformat_header($self->{vcf});
}

sub _add_headers {
  my $self = shift;
  my $vcf = $self->{'vcf'};

  my %options = (
    input => basename($self->{input}),
    ref => basename($self->{ref}),
    # probably more if we expose cutoffs
  );
  Sanger::CGP::Vcf::VcfUtil::add_vcf_process_log($vcf,
    Sanger::CGP::Vcf::VcfProcessLog->new(
        -input_vcf_source => basename($0),
        -input_vcf_ver => Sanger::CGP::Pindel->VERSION,
        -input_vcf_param => \%options,
      )
  );

  $vcf->add_header_line({'key'=>'source', 'value' => basename($0)}, 'append' => 1);

  my @format = (
    {key => 'FORMAT', ID => 'WTP', Number => 1, Type => 'Integer', Description => q{+ve strand reads BLATed to reference sequence at this location, input alignment depth when WTM='.'}},
    {key => 'FORMAT', ID => 'WTN', Number => 1, Type => 'Integer', Description => q{-ve strand reads BLATed to reference sequence at this location, input alignment depth when WTM='.'}},
    {key => 'FORMAT', ID => 'WTM', Number => 1, Type => 'Float', Description => q{Mismatch fraction of reads BLATed to reference sequence at this location (3 d.p.), '.' when no reads found via BLAT}},
    {key => 'FORMAT', ID => 'MTP', Number => 1, Type => 'Integer', Description => q{+ve strand reads BLATed to alternate sequence at this location}},
    {key => 'FORMAT', ID => 'MTN', Number => 1, Type => 'Integer', Description => q{-ve strand reads BLATed to alternate sequence at this location}},
    {key => 'FORMAT', ID => 'MTM', Number => 1, Type => 'Float', Description => q{Mismatch fraction of reads BLATed to alternate sequence at this location (3 d.p.), '.' when no reads found via BLAT}},
    {key => 'FORMAT', ID => 'VAF', Number => 1, Type => 'Float', Description => q{Variant allele fraction using reads that unambiguously map to ref or alt seq (3 d.p.)'}},
  );
  $self->{fmt_ext} = q{};
  for my $f(@format) {
    $vcf->add_header_line($f);
    $self->{fmt_ext} .= q{:}.$f->{ID};
  }
}

sub to_data_hash {
  my ($self, $v_d) = @_;
  my @items = @{$v_d};
  my %out;
  # simplified version of Vcf->next_data_hash
  # only stuff we need
  $out{CHROM}  = $items[0];
  $out{POS}    = $items[1];
  # trim the first base from these
  $out{REF}    = substr $items[3], 1;
  $out{ALT}    = substr $items[4], 1;

  # parse the info block
  for my $info (split(/;/,$items[7])) {
    my ($key,$val) = split(/=/,$info);
    # all the values we need are key/val
    next unless(defined $val);
    die "Clash between INFO and columns '$key'" if(exists $out{$key});
    $out{$key} = $val;
  }

  # add END
  $out{END} = $out{POS};
  if($out{PC} ne 'I') { # so D/DI
    $out{END} += $out{LEN} + 1;
  }
  else {
    $out{END} += 1;
  }

  # skip FORMAT
  # parse GT
  my ($gt, $pp, $pn) = split /:/, $items[9];
  # put in top level as no clash
  $out{PP} = $pp;
  $out{PN} = $pn;
  return \%out;
}

sub process_records {
  my $self = shift;
  my $fh = $self->{ofh};
  my $readtmp_dir = tempdir( CLEANUP => 0 ); # doesn't clear as you would expect
  my $last_chr_pos = q{.};
  while(my $v_d = $self->{vcf}->next_data_array) {
    my $this_c_p = sprintf "%s:%d", $v_d->[0], $v_d->[1];
    if($last_chr_pos ne $this_c_p) {
      # we use very large search range, so if the start pos doesn't change don't want to reparse reads
      remove_tree($readtmp_dir, { keep_root => 1 });
      $last_chr_pos = $this_c_p;
    }
    $v_d->[$V_FMT] .= $self->{fmt_ext} unless($self->{fill_in});
    $self->blat_record($v_d, $readtmp_dir);
    printf $fh "%s\n", join "\t", @{$v_d};
  }
  if(-d $readtmp_dir) {
    remove_tree($readtmp_dir);
  }
  $self->_close_sams
}

sub blat_record {
  my ($self, $v_d, $readtmp_dir) = @_;
  my $v_h = $self->to_data_hash($v_d);

  # now attempt the blat stuff
  my ($fh_target, $file_target) = tempfile( SUFFIX => '.fa', UNLINK => 1 );
  $self->blat_ref_alt($fh_target, $v_h);
  close $fh_target or die "Failed to close blat ref temp file";

  my $change_pos_low = $v_h->{change_pos};
  $change_pos_low++ if($v_h->{PC} eq 'I');
  my $range_l = ($v_h->{RE} - $v_h->{RS}) + 1;
  my $change_pos_high = $change_pos_low + $range_l; # REF based range, adjusted in func
  $v_h->{change_pos_low} = $change_pos_low;
  $v_h->{change_pos_high} = $change_pos_high;

  my $gt_pos = $V_GT_START-1;
  for my $sample(@{$self->{vcf_sample_order}}) {
    $gt_pos++;
    if($self->{fill_in} && $v_d->[$gt_pos] ne q{.}) {
      next;
    }
    my $gt_set = $self->blat_reads($v_h, $file_target, $readtmp_dir, $sample);
    if($v_d->[$gt_pos] eq q{.}) {
      $v_d->[$gt_pos] = join q{:}, './.:.:.:.:.', @{$gt_set};
    }
    else {
      $v_d->[$gt_pos] = join q{:}, $v_d->[$gt_pos], @{$gt_set};
    }
  }
  # tempfile unlink only does it on shutdown when used in this way
  unlink $file_target;
  return 1;
}

sub read_ranges {
  my ($self, $v_h, $sample) = @_;
  # return a string of chr:s-e... if approprate.
  my $read_buffer = $self->{max_insert};
  return sprintf $LOCI_FMT, $v_h->{CHROM}, $v_h->{q_start} - $read_buffer, $v_h->{q_end} + $read_buffer;
}

sub blat_reads {
  my ($self, $v_h, $file_target, $readtmp_dir, $sample) = @_;
  my $file_query = sprintf '%s/%s.fa', $readtmp_dir, $sample;
  # setup the temp file
  my ($fh_psl, $file_psl) = tempfile( SUFFIX => '.psl', UNLINK => 1);
  close $fh_psl or die "Failed to close $file_psl (psl output)";

  my $c_reads = sprintf $READS_ONLY, $self->{hts}->{$sample}->hts_path, $self->read_ranges($v_h, $sample), $file_query;
  if(! -e $file_query) {
    my ($r_out, $r_err, $r_exit) = capture { system([0], $c_reads); };
  }

  my $c_blat = sprintf $BLAT_ONLY, $file_target, $file_query, $file_psl, $file_psl, $file_target, $file_query;
  my ($c_out, $c_err, $c_exit) = capture { system([0,255], $c_blat); };
  if($c_exit == 255 && ($c_err =~ m/processed 0 reads/ms || $c_err =~ m/End of file reading 4 bytes/ms)) {
    # No reads found
    $c_exit = 0;
  }
  if($c_exit) {
    warn "An error occurred while executing: $c_blat\n";
    warn "\tERROR: $c_err\n";
    warn "\tECODE: $c_exit\n";
    warn "Read command: $c_reads\n";
    warn "DATA BLOCK\n";
    warn "Target:\n";
    my ($t_out, $t_err, $t_exit) = capture { system("cat $file_target"); };
    warn $t_out;
    warn "Query:\n";
    ($t_out, $t_err, $t_exit) = capture { system("cat $file_query"); };
    warn $t_out;
    warn "PSL:\n";
    ($t_out, $t_err, $t_exit) = capture { system("cat $file_psl"); };
    warn $t_out;
    exit $c_exit;
  }

  # tempfile unlink only does it on shutdown when used in this way
  unlink $file_psl;

  my ($wtp, $wtn, $mtp, $mtn, $wt_bmm, $mt_bmm) = $self->psl_axt_parser(\$c_out, $v_h, $sample);
  my ($wtm, $mtm) = (q{.}, q{.});
  my $wtr = $wtp + $wtn;
  if($wtr > 0) {
    $wtm = sprintf("%.3f", $wt_bmm / $wtr);
  }
  my $mtr = $mtp+$mtn;
  if($mtr > 0) {
    $mtm = sprintf("%.3f", $mt_bmm / $mtr);
  }
  my $depth = $wtr+$mtr;
  my $vaf = sprintf("%.3f", $depth ? $mtr/$depth : 0);
  return [$wtp, $wtn, $wtm, $mtp, $mtn, $mtm, $vaf];
}

sub psl_axt_parser {
  my ($self, $blat_axt, $v_h, $sample) = @_;
  # collate the data by readname and order by score
  my @lines = split /\n/, ${$blat_axt};
  my $line_c = @lines;
  # group all reads and order by score
  my %reads;
  for(my $i = 0; $i<$line_c; $i+=4) {
    my ($id, $t_name, $t_start, $t_end, $q_name, $q_start, $q_end, $strand, $score) = split q{ }, $lines[$i];
    next if($score < 0); # it happens
    my $clean_qname = $q_name;
    $clean_qname =~ s{/([12])$}{};
    my $q_seq = $lines[$i+2];
    next if(length $q_seq < int $self->{min_rl} * $MAPPED_RL_MULT);
    push @{$reads{$clean_qname}{$score}}, [$t_name, $t_start, $t_end, $q_name, $q_start, $q_end, $strand, $score, $lines[$i+1], $q_seq];
  }

  my %type_strand;
  my %bmm_sums;
  my @ref_reads;
  my $is_del = 0;
  if(length $v_h->{change_ref} > length $v_h->{change_alt}) {
    $is_del = 1;
  }
  # sort keys for consistency
  READ: for my $read(sort keys %reads) {
    # get the alignment with the highest score
    for my $score(sort {$b<=>$a} keys %{$reads{$read}}) {
      my @records = @{$reads{$read}{$score}};
      if(@records != 1) { # if best score has more than one alignment it is irrelevant
        next READ;
      }
      my $record = $records[0];
      my $ref_or_alt = $record->[0];
      if($is_del == 1 && $ref_or_alt eq 'REF') {
        # see block at end of this function
        push @ref_reads, $record;
      }
      if($self->parse_axt_event($v_h, $record, $sample) == 1) {
        $type_strand{$ref_or_alt.$record->[6]} += 1;
        $bmm_sums{$ref_or_alt} += bmm($record->[8], $record->[9]);
      }
      next READ; # remaining items are worse alignments
    }
  }
  my $wtp = $type_strand{'REF+'} || 0;
  my $wtn = $type_strand{'REF-'} || 0;
  my $mtp = $type_strand{'ALT+'} || 0;
  my $mtn = $type_strand{'ALT-'} || 0;

  # only relevant for deletions
  if($is_del == 1 && $wtp == 0 && $wtn == 0 && $v_h->{RE} - $v_h->{RS} > $LARGE_D) {
    # if a deletion is large it can cause no REF depth as impossible for a read to span the ends of the event.
    # rather than specifying a cutoff we rely on the data to drive this
    my $add_bmb_sum;
    ($wtp, $wtn, $add_bmb_sum) = $self->parse_axt_del_ref($v_h, \@ref_reads, $sample);
    $bmm_sums{'REF'} += $add_bmb_sum;
  }

  return ($wtp, $wtn, $mtp, $mtn, $bmm_sums{REF}, $bmm_sums{ALT});
}

sub parse_axt_del_ref {
  my ($self, $v_h, $records, $sample) = @_;
  # Need to return
  #   pos reads
  #   neg reads
  #   mismatch fractions (via bmm)
  # Augment $self->sam_record($v_h, $rec, $sample);

  # some rules:
  # - this shouldn't be getting called if the event spans both ends so can assume reads at each end can be saved without a check
  # - require near perfect match to ref, but reads that get here are known NOT to have a better ALT mapping

  my ($wtp, $wtn, $bmm_sum) = (0,0,0);
  my $change_len = $v_h->{change_pos_high} - $v_h->{change_pos_low};
  my $mid_point = $v_h->{change_pos_low} + int($change_len/2);
  for my $rec(@{ $records }) {
    my ($t_name, $t_start, $t_end, $q_name, $q_start, $q_end, $strand, $score, $t_seq, $q_seq) = @{ $rec };
    if(
        ($t_start < $v_h->{change_pos_low} && $t_end > $v_h->{change_pos_low})
        ||
        ($t_start < $v_h->{change_pos_high} && $t_end > $v_h->{change_pos_high})
        #($t_start < $v_h->{change_pos_low} && $t_end > $v_h->{change_pos_low})
    ) {
      if($strand eq '+') {
        $wtp += 1;
      }
      else {
        $wtn += 1;
      }
      $bmm_sum += bmm($t_seq, $q_seq);
      $self->sam_record($v_h, $rec, $sample);
    }
  }
  # we've counted 2 locations so can't just return the full value
  # this is under review
  if($wtp > 0) {
    $wtp = int($wtp/2) + 1;
  }
  if($wtn > 0) {
    $wtn = int($wtn/2) + 1;
  }
  return ($wtp, $wtn, $bmm_sum);
}

sub parse_axt_event {
  my ($self, $v_h, $rec, $sample) = @_;
  my ($t_name, $t_start, $t_end, $q_name, $q_start, $q_end, $strand, $score, $t_seq, $q_seq) = @{$rec};

  # specific to deletion class
  my $change_seq = $v_h->{change_ref};
  my $change_pos_high = $v_h->{change_pos_high};
  if($t_name eq 'ALT') {
    $change_pos_high -= $v_h->{LEN};
    $change_seq = $v_h->{change_alt};
  }

  # all the reads that span the range are kept
  my $retval = 0;
  if($t_start <= ($v_h->{change_pos_low} - $PAD_EVENT) && $t_end > ($change_pos_high + $PAD_EVENT)) {
    # look for the change (or absence) where we expect it
    my $exp_pos = $v_h->{change_pos_low} - $t_start;
    if($exp_pos <= length $q_seq) {
      my $sub_q_seq = substr($q_seq, $exp_pos, length $change_seq);
      if(length $change_seq == length $sub_q_seq # same length
        && index($t_seq, q{-}) == -1 # no gaps
        && index($q_seq, q{-}) == -1 # no gaps
        && substr($change_seq,0,1) eq substr($sub_q_seq,0,1) # matching first base
        && substr($change_seq,-1,1) eq substr($sub_q_seq,-1,1) # matching last base
      ) {
        $retval = 1;
        $self->sam_record($v_h, $rec, $sample);
      }
    }
  }
  return $retval;
}

sub bmm {
  my ($a, $b) = @_; # need a copy of the strings anyway
  my $len = length $a;
  my $diffs = 0;
  for(0..($len-1)) {
    $diffs++ if(chop $a ne chop $b);
  }
  return $diffs/$len;
}

sub sam_record {
  my($self, $v_h, $rec, $sample) = @_;

  my $qname = $rec->[3];
  my $seq = $rec->[9];
  my $flag = 0; # not paired
  $flag += 16 if($rec->[6] eq '-');

  # POS is the base preceeing any change, seq start it this - target_pad
  my $pos = ($v_h->{POS} - $self->{target_pad}) + $rec->[1];
  my $cigar = q{};
  if($rec->[0] eq 'REF') {
    $cigar = length($seq).'M';
  }
  else {
    my $m_c = ($v_h->{change_pos} - $rec->[1]) + 1;
    $m_c += 1 if($v_h->{PC} eq 'I');
    $cigar = $m_c.'M';
    my $change_ref = length($v_h->{REF});
    my $change_alt = length($v_h->{ALT});
    if($change_ref) {
      $cigar .= $change_ref.'D';
    }
    if($change_alt) {
      $cigar .= $change_alt.'I';
      $m_c += $change_alt; # as consumes read
    }
    $cigar .= (length($seq) - $m_c).'M';
  }
  printf {$self->{sfh}->{$sample}} "%s\n", join "\t", $qname, $flag, $v_h->{CHROM}, $pos, 60, $cigar, '*', 0, 0, $seq, '*';
}

sub sam_to_bam {
  my ($self) = @_;
  for my $sample(@{$self->{vcf_sample_order}}) {
    my $sam = $self->{samfile}->{$sample};
    my $bam = $self->{bamfile}->{$sample};
    my $tmp = $bam;
    $tmp =~ s/bam$/tmp/;
    my $command = sprintf q{bash -c 'set -o pipefail ; zcat %s | pee "grep ^@" "grep -v ^@ | sort | uniq" | samtools view -uT %s - | samtools sort -l 0 -T %s - | samtools calmd - %s > %s'},
                  $sam, # zcat
                  $self->{ref}, # view
                  $tmp, # sort
                  $self->{'ref'}, $bam; # calmd
    my ($c_out, $c_err, $c_exit) = capture { system($command); };
    if($c_exit) {
      warn "An error occurred while executing $command\n";
      warn "\tERROR$c_err\n";
      exit $c_exit;
    }
    unlink(glob(sprintf '%s.*.bam', $tmp));
    unlink $sam;
  }
}

sub flanking_ref {
  my ($self, $v_h) = @_;
  my $ref_left = $self->{fai}->get_sequence_no_length(
    sprintf $LOCI_FMT,
    $v_h->{CHROM},
    ($v_h->{POS} - $self->{target_pad})+1,
    $v_h->{POS},
  );

  my $ref_right = $self->{fai}->get_sequence_no_length(
    sprintf $LOCI_FMT,
    $v_h->{CHROM},
    $v_h->{END},
    $v_h->{END} + $self->{target_pad},
  );

  return [$ref_left, $ref_right]
}

sub blat_ref_alt {
  my ($self, $fh, $v_h) = @_;
  my ($ref_left, $ref_right) = @{$self->flanking_ref($v_h)};
  my $ref = $v_h->{REF};
  my $alt = $v_h->{ALT};
  my $r_start = $v_h->{RS};
  my $r_end = $v_h->{RE};
  my $change_at = length $ref_left;

  my $call_type = $v_h->{PC};
  if($call_type eq 'I') {
    $change_at -= 1; # force base before
  }

  # used when getting reads from HTSfile
  my $q_start = $v_h->{POS} - $change_at; # correcting for position handled in change_at
  my $q_end = $v_h->{POS} + length $ref_right;
  if($call_type eq 'D') {
    $q_end += length $v_h->{REF};
  }

  print $fh sprintf ">REF\n%s%s%s\n", $ref_left, $ref, $ref_right or die "Failed to write REF to blat ref temp file";
  print $fh sprintf ">ALT\n%s%s%s\n", $ref_left, $alt, $ref_right or die "Failed to write ALT to blat ref temp file";

  my $seq_left = substr($ref_left, -1);

  # -1 as includes the base before and after which would be -2 but need to correct for coord maths
  # (for Del and Ins, unsure about DI at the moment)
  my $seq_right;
  if($call_type ne 'I') {
    $seq_right = substr($ref_right, 0, ($r_end - $r_start) - $v_h->{LEN});
  }
  else {
    $seq_right = substr($ref_right, 0, ($r_end - $r_start));
  }

  my $change_ref = $seq_left.$ref.$seq_right;
  my $change_alt = $seq_left.$alt.$seq_right;

  $v_h->{q_start} = $q_start;
  $v_h->{q_end} = $q_end;
  $v_h->{change_pos} = $change_at;
  $v_h->{change_ref} = $change_ref;
  $v_h->{change_alt} = $change_alt;
  return 1;
}
