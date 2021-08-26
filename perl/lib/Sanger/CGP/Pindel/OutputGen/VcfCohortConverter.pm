# Copyright (c) 2014-2021
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
package Sanger::CGP::Pindel::OutputGen::VcfCohortConverter;
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
    {key => 'INFO', ID => 'GC5P', Number => 1, Type => 'Float', Description => 'GC content of 200 bp. 5 prime'},
    {key => 'INFO', ID => 'GCRNG', Number => 1, Type => 'Float', Description => 'GC content of deleted/inserted seq, including range'},
    {key => 'INFO', ID => 'GC3P', Number => 1, Type => 'Float', Description => 'GC content of 200 bp. 3 prime'}
  );

  my @format = (
    {key => 'FORMAT', ID => 'GT', Number => 1, Type => 'String', Description => 'Genotype'},
    {key => 'FORMAT', ID => 'S1', Number => 1, Type => 'Integer', Description => 'Pindel S1 score'},
    {key => 'FORMAT', ID => 'S2', Number => 1, Type => 'Float', Description => 'Pindel S2 score, not present for all types'},
    {key => 'FORMAT', ID => 'PP', Number => 1, Type => 'Integer', Description => 'Pindel calls on the positive strand'},
    {key => 'FORMAT', ID => 'NP', Number => 1, Type => 'Integer', Description => 'Pindel calls on the negative strand'},
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
  $ret .= 'PC='.$record->type().';';
  $ret .= 'RS='.$record->range_start().';';
  $ret .= 'RE='.$record->range_end().';';
  $ret .= 'LEN='.$record->length().';';
  $ret .= 'REP='.$record->repeats().';';
  $ret .= sprintf 'GC5P=%.3f;', $record->gc_5p;
  $ret .= sprintf 'GCRNG=%.3f;', $record->gc_rng;
  $ret .= sprintf 'GC3P=%.3f', $record->gc_3p;
  $ret .= $SEP;

  # FORMAT
  $ret .= $self->{_format};

  for my $samp(@{$self->{_srt_samples}}) {
    $ret .= $SEP;
    if($self->gen_all || exists $record->reads->{$samp}) {
      $ret .= './.:';
      $ret .= $record->s1.q{:};
      $ret .= $record->s2 || '.';
      $ret .= q{:};
      $ret .= $record->get_read_counts($samp, '+').q{:};
      $ret .= $record->get_read_counts($samp, '-');
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

sub _dump_rec_detail {
  my $self = shift;
  for my $k(sort qw(change_pos_low change_pos_high change_l change_ref change_alt q_start q_end type chr)) {
    printf "%s: %s\n", $k, ref $self->{$k} ? Dumper($self->{$k}) : $self->{$k};
  }
  print "\n";
  return 0;
}
