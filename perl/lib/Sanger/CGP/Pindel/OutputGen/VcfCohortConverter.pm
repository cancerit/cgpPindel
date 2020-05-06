package Sanger::CGP::Pindel::OutputGen::VcfCohortConverter;

########## LICENCE ##########
# Copyright (c) 2020 Genome Research Ltd.
#
# Author: CASM IT <cgphelp@sanger.ac.uk>
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
use File::Basename;
use File::Temp qw(tempfile);
use Capture::Tiny qw(capture);

use Data::Dumper;

use Sanger::CGP::Pindel;
use Sanger::CGP::Vcf::VcfUtil;
use Sanger::CGP::Vcf::VcfProcessLog;
use Const::Fast qw(const);
use Sanger::CGP::Pindel::InputGen;

use Vcf;

const my $SEP => "\t";
const my $NL => "\n";
const my $READS_AND_BLAT => q{bash -c 'set -o pipefail ; samtools view -uF 3840 %s %s:%d-%d | samtools fasta - > %s && blat -t=dna -q=dna -noTrimA -minIdentity=95 -noHead -out=psl %s %s %s && pslPretty -axt %s %s %s /dev/stdout'};
# returns +ve and then -ve results
const my $SAM_DEPTH_PN => q{bash -c "set -o pipefail ; samtools view -uF 3844 %s %s:%d-%d | pee 'samtools view -c -F 16 -' 'samtools view -c -f 16 -'"};

const my $MATCH => 0;
const my $Q_GAP => 4;
const my $T_GAP => 6;
const my $T_BASES => 7;
const my $STRAND => 8;
const my $Q_NAME => 9;
const my $Q_SIZE => 10;
const my $Q_START => 11;
const my $T_NAME => 13;
const my $T_START => 15;
const my $T_END => 16;
const my $BLOCK_COUNT => 17;
const my $BLOCK_SIZES => 18;
const my $Q_STARTS => 19;
const my $T_STARTS => 20;
const my $Q_SEQ => 21;
const my $T_SEQ => 22;

const my $SD_MULT => 2;

1;

sub new{
  my $proto = shift;
  my (%args) = @_;
  my $class = ref($proto) || $proto;

  my $self = {};
  bless $self, $class;

  $self->init(%args);

  return $self;
}

sub init{
  my($self,%args) = @_;
  $self->{_contigs} = $args{-contigs};
  $self->{_srt_samples} = $args{-samples};
  $self->{_hts_set} = $args{-hts_set};
  $self->{_bas_set} = $args{-bas_set};
  $self->_max_inserts() if(defined $self->{_bas_set});
  $self->{_all} = $args{-all};
  if(defined $args{-badloci}) {
    $self->{_tabix} = Sanger::CGP::Pindel::InputGen::tabix_to_interval_tree($args{-badloci});
  }
}

sub _interval_hit {
  my ($self, $chr, $start, $stop) = @_;
  return 0 unless (exists $self->{_tabix}->{$chr});
  return scalar @{$self->{_tabix}->{$chr}->fetch($start, $stop)};
}

sub _max_inserts {
  my $self = shift;
  my %bas = %{$self->{_bas_set}};
  my %ins_by_sample;
  for my $s(keys %bas) {
    my $max_ins = 0;
    for my $rg($bas{$s}->read_groups) {
      my $m_sd = $bas{$s}->get($rg, 'mean_insert_size') + ($bas{$s}->get($rg, 'insert_size_sd') * $SD_MULT);
      $max_ins = $m_sd if($m_sd > $max_ins);
    }
    $ins_by_sample{$s} = $max_ins;
  }
  $self->{_ins_set} = \%ins_by_sample;
}


=head gen_header

Generates a Vcf header String for NORMAL/TUMOUR comparisons.

@param1 reference_path - a String containing the path to the reference used in the VCF.

@param2 input_source   - a String containing the name and version of the application or source of the VCF data.

@param3 sample         - hash-ref of a Sanger::CGP::Vcf::Sample objects representing samples to be included.

@param3 options        - hash-ref of options passed to generating command

=cut
sub gen_header{
  my($self, $reference_path, $input_source, $samples, $options) = @_;

  my @process_logs = (
    Sanger::CGP::Vcf::VcfProcessLog->new(
      -input_vcf_source => 'Pindel',
      -input_vcf_ver => 'v02', # always have S2 at this point
    ),
    Sanger::CGP::Vcf::VcfProcessLog->new(-input_vcf_source => basename($0),
      -input_vcf_ver => Sanger::CGP::Pindel->VERSION,
      -input_vcf_param => $options,
    ),
  );

  my @info = (
    {key => 'INFO', ID => 'PC', Number => 1, Type => 'String', Description => 'Pindel call'},
    {key => 'INFO', ID => 'RS', Number => 1, Type => 'Integer', Description => 'Range start'},
    {key => 'INFO', ID => 'RE', Number => 1, Type => 'Integer', Description => 'Range end'},
    {key => 'INFO', ID => 'LEN', Number => 1, Type => 'Integer', Description => 'Length'},
    {key => 'INFO', ID => 'REP', Number => 1, Type => 'Integer', Description => 'Change repeat count within range'},
    {key => 'INFO', ID => 'S1', Number => 1, Type => 'Integer', Description => 'S1'},
    {key => 'INFO', ID => 'S2', Number => 1, Type => 'Float', Description => 'S2'}
  );

  my @format = (
    {key => 'FORMAT', ID => 'GT', Number => 1, Type => 'String', Description => 'Genotype'},
    {key => 'FORMAT', ID => 'S1', Number => 1, Type => 'Integer', Description => 'Pindel S1 score'},
    {key => 'FORMAT', ID => 'S2', Number => 1, Type => 'Float', Description => 'Pindel S2 score, not present for all types'},
    {key => 'FORMAT', ID => 'PP', Number => 1, Type => 'Integer', Description => 'Pindel calls on the positive strand'},
    {key => 'FORMAT', ID => 'NP', Number => 1, Type => 'Integer', Description => 'Pindel calls on the negative strand'},
    #{key => 'FORMAT', ID => 'WTP', Number => 1, Type => 'Integer', Description => '+ve strand reads BLATed (or count for large deletions) to reference at this location'},
    #{key => 'FORMAT', ID => 'WTN', Number => 1, Type => 'Integer', Description => '-ve strand reads BLATed (or count for large deletions) to reference at this location'},
    #{key => 'FORMAT', ID => 'MTP', Number => 1, Type => 'Integer', Description => '+ve strand reads BLATed to alternate sequence at this location'},
    #{key => 'FORMAT', ID => 'MTN', Number => 1, Type => 'Integer', Description => '-ve strand reads BLATed to alternate sequence at this location'},
    #{key => 'FORMAT', ID => 'VAF', Number => 1, Type => 'Float', Description => 'Variant allele fraction using reads that unabiguously map to ref or alt seq (to 3 d.p.)'},
  );

  my @blank_fmt = (q{.}) x (scalar @format -1);
  $self->{_noread_gt} = join q{:}, './.', @blank_fmt;

  my $fmt_str = q{};
  for my $f(@format) {
    $fmt_str .= q{:} if($fmt_str);
    $fmt_str .= $f->{ID};
  }
  $self->{_format} = $fmt_str;

  my $vcf = Vcf->new(version=>'4.1');
  my @timeData = localtime(time);
  $vcf->add_header_line( { key => 'fileDate', value => sprintf '%d%02d%02d', 900 + $timeData[5], $timeData[4]+1, $timeData[3] } );
  $vcf->add_header_line( { key => 'source',   value => $input_source }, 'append' => 1 );
  $vcf->add_header_line( { key => 'reference', value => $reference_path } );

  for my $contig (@{$self->{_contigs}}){
    Sanger::CGP::Vcf::VcfUtil::add_vcf_contig($vcf,$contig)
  }

  for my $inf (@info){
    $vcf->add_header_line($inf);
  }

  for my $for (@format){
    $vcf->add_header_line($for);
  }

  for my $process_log (@process_logs){
    Sanger::CGP::Vcf::VcfUtil::add_vcf_process_log($vcf,$process_log)
  }

  for my $samp(sort keys %{$samples}) {
    Sanger::CGP::Vcf::VcfUtil::add_vcf_sample($vcf, $samples->{$samp}, $samp);
  }

  return $vcf->format_header();
}

sub gen_record{
  my($self, $record) = @_;
  my $ret = q{};

  # CHR POS ID REF ALT QUAL FILTER INFO FORMAT GENO GENO

  my $start = $record->start();
  $start-- if(substr($record->type(),0,1) eq 'D');

  my $ref = uc ($record->lub . $record->ref_seq);
  my $alt = uc ($record->lub . $record->alt_seq);

  if($self->_interval_hit($record->chro, $start, $start + length $ref)) {
    return $ret;
  }

  $ret .= $record->chro().$SEP;
  $ret .= $start.$SEP;
  $ret .= $record->id().$SEP;

  $ret .= $ref.$SEP;
  $ret .= $alt.$SEP;
  $ret .= $record->sum_ms().$SEP;
  $ret .= '.'.$SEP;

  # INFO
  #PC=D;RS=19432;RE=19439;LEN=3;S1=4;S2=161.407;REP=2;PRV=1
  $ret .= 'PC='.$record->type().';';
  $ret .= 'RS='.$record->range_start().';';
  $ret .= 'RE='.$record->range_end().';';
  $ret .= 'LEN='.$record->length().';';
  $ret .= 'REP='.$record->repeats().$SEP;

  # # now attempt the blat stuff
  # my ($fh_target, $file_target) = tempfile(
  #   DIR => '/dev/shm',
  #   SUFFIX => '.fa',
  #   UNLINK => 1
  # );

  # my ($q_start, $q_end, $change_pos, $change_ref, $change_alt) = $self->blat_ref_alt($fh_target, $record);

  # my $change_pos_low = $change_pos;
  # $change_pos_low++ if($record->type eq 'I');
  # my $range_l = ($record->range_end - $record->range_start) + 1;
  # my $change_pos_high = $change_pos_low + $range_l; # REF based range, adjusted in func
  # my $change_l = $record->end - $record->start + 1;
  # $self->{change_pos_low} = $change_pos_low;
  # $self->{change_pos_high} = $change_pos_high;
  # $self->{change_l} = $change_l;
  # $self->{change_ref} = $change_ref;
  # $self->{change_alt} = $change_alt;
  # $self->{q_start} = $q_start;
  # $self->{q_end} = $q_end;
  # $self->{file_target} = $file_target;
  # $self->{type} = $record->type;
  # $self->{chr} = $record->chro;

  # FORMAT
  $ret .= $self->{_format};

  for my $samp(@{$self->{_srt_samples}}) {
    #my ($wtp, $wtn, $mtp, $mtn) = $self->blat_reads($samp, $record);
    $ret .= $SEP;
    if($self->gen_all || exists $record->reads->{$samp}) {
      $ret .= './.:';
      $ret .= $record->s1.q{:};
      $ret .= $record->s2 || '.';
      $ret .= q{:};
      $ret .= $record->get_read_counts($samp, '+').q{:};
      $ret .= $record->get_read_counts($samp, '-');
      #.q{:};
      #$ret .= $wtp.q{:};
      #$ret .= $wtn.q{:};
      #$ret .= $mtp.q{:};
      #$ret .= $mtn.q{:};
      #my $mtr = $mtp+$mtn;
      #my $depth = $wtp+$wtn+$mtr;
      #$ret .= $depth ? sprintf("%.3f", $mtr/$depth) : 0;
    }
    else {
      $ret .= $self->{_noread_gt};
    }
  }
  $ret .= $NL;
  return $ret;
}

sub gen_all {
  return shift->{_all};
}

sub hts_file_by_sample {
  my ($self, $sample) = @_;
  return $self->{_hts_set}->{$sample};
}

sub bas_by_sample {
  my ($self, $sample) = @_;
  return $self->{_bas_set}->{$sample};
}

sub ins_by_sample {
  my ($self, $sample) = @_;
  return $self->{_ins_set}->{$sample};
}

sub TO_REMOVE_blat_reads {
  my ($self, $sample, $record) = @_;
  my ($fh_query, $file_query) = tempfile(
    DIR => '/dev/shm',
    SUFFIX => '.fa',
    UNLINK => 1
  );
  close $fh_query or die "Failed to close $file_query (query reads)";
  my ($fh_psl, $file_psl) = tempfile(
    DIR => '/dev/shm',
    SUFFIX => '.psl',
    UNLINK => 1
  );
  close $fh_psl or die "Failed to close $file_psl (psl output)";


  my $read_buffer = $self->ins_by_sample($sample);
  my $c_blat = sprintf $READS_AND_BLAT, $self->hts_file_by_sample($sample), $self->{chr}, $self->{q_start}-$read_buffer, $self->{q_end}+$read_buffer, $file_query, $self->{file_target}, $file_query, $file_psl, $file_psl, $self->{file_target}, $file_query;
  my ($c_out, $c_err, $c_exit) = capture { system($c_blat); };
  if($c_exit) {
    warn "An error occurred while executing $c_blat\n";
    warn "\tERROR$c_err\n";
    exit $c_exit;
  }

# print "target.fa\n";
# system("cat $self->{file_target}");
# <STDIN>;
# print "query.fa\n";
# system("cat $file_query");
# <STDIN>;
# print "\n$c_out\n";
# exit 1;

  #my ($wtp, $wtn, $mtp, $mtn) = $self->blat_counts(\$c_out, $sample);
  my ($wtp, $wtn, $mtp, $mtn) = $self->psl_axt_parser(\$c_out, $sample);
  if($wtp + $wtn == 0) {
    ($wtp, $wtn) = $self->sam_depth($sample, $record);
  }
  return ($wtp, $wtn, $mtp, $mtn);
}

sub TO_REMOVE_psl_axt_parser {
  my ($self, $blat_axt) = @_;
  # collate the data by readname and order by score
  my @lines = split /\n/, ${$blat_axt};
  my $line_c = @lines;
  # group all reads and order by score
  my %reads;
  for(my $i = 0; $i<$line_c; $i+=4) {
    my ($id, $t_name, $t_start, $t_end, $q_name, $q_start, $q_end, $strand, $score) = split q{ }, $lines[$i];
    push @{$reads{$q_name}{$score}}, [$t_name, $t_start, $t_end, $q_name, $q_start, $q_end, $strand, $score, $lines[$i+1], $lines[$i+2]];
  }
  my %TYPE_STRAND;
  # sort keys for consistency
  READ: for my $read(sort keys %reads) {
    # get the alignment with the highest score
    for my $score(sort {$b<=>$a} keys %{$reads{$read}}) {
      my @records = @{$reads{$read}{$score}};
      next READ if(@records != 1); # if best score has more than one alignment it is irrelevant
      my $record = $records[0];
      if($self->parse_axt_event($record) == 1) {
        $TYPE_STRAND{$record->[0].$record->[6]} += 1;
      }
      next READ; # remaining items are worse alignments
    }
  }
  my $wtp = $TYPE_STRAND{'REF+'} || 0;
  my $wtn = $TYPE_STRAND{'REF-'} || 0;
  my $mtp = $TYPE_STRAND{'ALT+'} || 0;
  my $mtn = $TYPE_STRAND{'ALT-'} || 0;
  return ($wtp, $wtn, $mtp, $mtn);
}

sub TO_REMOVE_parse_axt_event {
  my ($self, $rec) = @_;
  my ($t_name, $t_start, $t_end, $q_name, $q_start, $q_end, $strand, $score, $t_seq, $q_seq) = @{$rec};

  # specific to deletion class
  my $change_seq = $self->{change_ref};
  my $change_pos_high = $self->{change_pos_high};
  if($t_name eq 'ALT') {
    $change_pos_high -= $self->{change_l};
    $change_seq = $self->{change_alt};
  }
  $self->{change_seq} = $change_seq;

  # all the reads that don't span the range
  return 0 unless($t_start <= $self->{change_pos_low} && $t_end > $change_pos_high);
  # look for the change (or absence) where we expect it
  my $exp_pos = $self->{change_pos_low} - $t_start;
  if(index($t_seq, $change_seq, $exp_pos) == $exp_pos) {
# print join "\t", $t_name, $t_start, $t_end, $q_name, $q_start, $q_end, $strand, $score;
# printf "\n%s\n%s\n", $t_seq, $q_seq;
# print substr($t_seq, ($self->{change_pos_low} - $t_start), $change_pos_high)."\n";
# print "$change_seq\n";
# print "EXP: $exp_pos\n";
# print "INDEX: ".index($t_seq, $change_seq, $exp_pos)."\n";
    return 1;
  }

#<STDIN>;
#  return 2;
  return 0;
}


sub TO_REMOVE_sam_depth {
  my ($self, $sample, $record) = @_;
  my $mid_point = int ($record->range_start + (($record->range_end - $record->range_start)*0.5));
  my $c_samcount = sprintf $SAM_DEPTH_PN, $self->hts_file_by_sample($sample), $record->chro, $mid_point, $mid_point;
  my ($c_out, $c_err, $c_exit) = capture { system($c_samcount); };
  if($c_exit) {
    warn "An error occurred while executing $c_samcount\n";
    warn "\tERROR$c_err\n";
    exit $c_exit;
  }
  return (split /\n/, $c_out);
}

sub _dump_rec_detail {
  my $self = shift;
  for my $k(sort qw(change_pos_low change_pos_high change_l change_ref change_alt q_start q_end type chr)) {
    printf "%s: %s\n", $k, ref $self->{$k} ? Dumper($self->{$k}) : $self->{$k};
  }
  print "\n";
  return 0;
}

sub TO_REMOVE_blat_ref_alt {
  my ($self, $fh, $record) = @_;
  my $ref_left = $record->ref_left;
  my $ref_right = $record->ref_right;
  my $ref = q{};
  my $alt = q{}; # save an object lookup
  my $change_at = length $ref_left;

  if($record->type eq 'D') {
    $ref = uc $record->ref_seq;
    #nothing to add for alt
  } elsif($record->type eq 'I') {
    #nothing to add for ref
    $alt = $record->alt_seq;
    $change_at -= 1; # force base before
  }
  else { # must be DI
    $ref = $record->ref_seq;
    $alt = $record->alt_seq;
  }

  my $q_start = $record->start - $change_at; # correcting for position handled in change_at
  my $q_end = $record->end + length $ref_right;

  print $fh sprintf ">REF\n%s%s%s\n", $ref_left, $ref, $ref_right or die "Failed to write REF to blat ref temp file";
  print $fh sprintf ">ALT\n%s%s%s\n", $ref_left, $alt, $ref_right or die "Failed to write ALT to blat ref temp file";
  close $fh or die "Failed to close blat ref temp file";

  my $seq_left = substr($ref_left, -1);
  # -1 as includes the base before and after which would be -2 but need to correct for coord maths
  # (for Del and Ins, unsure about DI at the moment)
  my $seq_right;
  if($record->type eq 'D') {
    $seq_right = substr($ref_right, 0, ($record->range_end - $record->range_start) - 1);
  }
  elsif($record->type eq 'I') {
    $seq_right = substr($ref_right, 0, ($record->range_end - $record->range_start));
  }


  my $change_ref = uc ($seq_left.$ref.$seq_right);
  my $change_alt = uc ($seq_left.$alt.$seq_right);

# printf "%s\n", $record->ref_left;
# printf "%s\n", $record->ref_seq;
# printf "%s\n", $record->alt_seq;
# printf "%s\n", $record->ref_right;
# printf "%s\n", $record->type;
# printf "%s\n", $change_ref;
# printf "%s\n", $change_alt;
# #exit;

  return $q_start, $q_end, $change_at, $change_ref, $change_alt;
}
