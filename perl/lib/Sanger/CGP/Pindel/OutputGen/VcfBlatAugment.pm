package Sanger::CGP::Pindel::OutputGen::VcfBlatAugment;

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
use warnings FATAL => 'all';
use autodie qw(:all);
use Const::Fast qw(const);
use File::Basename;
use List::Util qw(max);
use File::Temp qw(tempfile);
use Capture::Tiny qw(capture);
use Vcf;
use Bio::DB::HTS;
use Bio::DB::HTS::Faidx;
use Sanger::CGP::Pindel;
use Sanger::CGP::Vcf::VcfProcessLog;
use Sanger::CGP::PindelPostProcessing::VcfSoftFlagger;
use Sanger::CGP::Vcf::VcfUtil;
use PCAP::Bam::Bas;

use Data::Dumper;

const my $SD_MULT => 2;
const my $V_FMT => 8;
const my $V_GT => 9;
const my $READS_AND_BLAT => q{bash -c 'set -o pipefail ; samtools view -uF 3840 %s %s | samtools fasta - > %s && blat -t=dna -q=dna -noTrimA -minIdentity=95 -noHead -out=psl %s %s %s && pslPretty -long -axt %s %s %s /dev/stdout'};
# returns +ve and then -ve results
const my $SAM_DEPTH_PN => q{bash -c "set -o pipefail ; samtools view -uF 3844 %s %s | pee 'samtools view -c -F 16 -' 'samtools view -c -f 16 -'"};
const my $LOCI_FMT => '%s:%d-%d';

1;

sub new{
  my $proto = shift;
  my (%args) = @_;
  my $class = ref($proto) || $proto;

  my $self = {
    input => $args{input},
    ref => $args{ref},
    ofh => $args{ofh},
    hts_file => $args{hts_file},
  };
  bless $self, $class;
  $self->_init;

  return $self;
}

sub _init {
  my $self = shift;
  my $vcf = Vcf->new(file => $self->{input});
  $vcf->parse_header();
  $self->{vcf} = $vcf;
  $self->_add_headers;
  $self->_hts;
  my $bas_sample = $self->_buffer_sizes; # has to go before validate sample
  $self->_validate_sample($bas_sample);
  # load the fai
  $self->{fai} = Bio::DB::HTS::Faidx->new($self->{ref});
  return 1;
}

sub _validate_sample {
  my ($self, $bas_sample) = @_;
  # check BAM/CRAM and VCF have same sample

  my $hts_sample;
  foreach my $line (split(/\n/,$self->{hts}->header->text)) {
    next unless($line =~ m/^\@RG/);
    chomp $line;
    ($hts_sample) = $line =~ m/SM:([^\t]+)/;
    last if(defined $hts_sample);
  }
  my @samples = $self->{vcf}->get_samples();
  die sprintf "ERROR: Only expecting 1 sample in VCF '%s', got %d\n", $self->{vcf}, scalar @samples if(@samples > 1);
  die sprintf "ERROR: Sample mismatch between BAM/CRAM (%s) and VCF (%s)\n", $hts_sample, $samples[0] if($hts_sample ne $samples[0]);
  die sprintf "ERROR: Sample mismatch between BAM/CRAM (%s) and BAS (%s)\n", $hts_sample, $bas_sample if($hts_sample ne $bas_sample);
  $self->{sample} = $hts_sample;
}

sub _buffer_sizes {
  my $self = shift;
  ## Set the max read length and insert size + SD*$SD_MULTI using the bas file data
  my $b = PCAP::Bam::Bas->new($self->{hts_file}.'.bas');
  my $max_ins = 0;
  my $max_rl = 0;
  my $sample;
  for my $rg($b->read_groups) {
    my $m_sd = $b->get($rg, 'mean_insert_size') + ($b->get($rg, 'insert_size_sd') * $SD_MULT);
    $max_ins = $m_sd if($m_sd > $max_ins);
    my $rl = max ($b->get($rg, 'read_length_r1'), $b->get($rg, 'read_length_r2'));
    $max_rl = $rl if($rl > $max_rl);
    my $s = $b->get($rg, 'sample');
    if($sample) {
      die "ERROR: Multiple samples found in %s\n", $self->{hts}.'.bas' if($sample ne $s)
    }
    else {
      $sample = $s;
    }
  }
  $self->{max_insert} = $max_ins;
  $self->{max_rl} = $max_rl;
  return $sample;
}

sub _hts {
  my $self = shift;
  $self->{hts} = Bio::DB::HTS->new(-bam => $self->{hts_file}, -fasta => $self->{ref});
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
    {key => 'FORMAT', ID => 'WTP', Number => 1, Type => 'Integer', Description => '+ve strand reads BLATed (or count for large deletions) to reference at this location'},
    {key => 'FORMAT', ID => 'WTN', Number => 1, Type => 'Integer', Description => '-ve strand reads BLATed (or count for large deletions) to reference at this location'},
    {key => 'FORMAT', ID => 'MTP', Number => 1, Type => 'Integer', Description => '+ve strand reads BLATed to alternate sequence at this location'},
    {key => 'FORMAT', ID => 'MTN', Number => 1, Type => 'Integer', Description => '-ve strand reads BLATed to alternate sequence at this location'},
    {key => 'FORMAT', ID => 'VAF', Number => 1, Type => 'Float', Description => 'Variant allele fraction using reads that unambiguously map to ref or alt seq (to 3 d.p.)'},
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
    die "Clash between INFO and columns" if(exists $out{$key});
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
#print "$out{POS} - $out{END}\n";
#exit;

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
my $count=0;
  while(my $v_d = $self->{vcf}->next_data_array) {
$count++;
#next if($count != 8); # important the one relating to the bug
#next if($count != 21247);
#printf "%s\n", join "\t", @{$v_d};
    my $v_h = $self->to_data_hash($v_d);
#next if($v_h->{'POS'} != 49092625);
    my @extra_gt = $self->blat_record($v_h);
    $v_d->[$V_FMT] .= $self->{fmt_ext};
    $v_d->[$V_GT] = join q{:}, $v_d->[$V_GT], @extra_gt;
    printf $fh "%s\n", join "\t", @{$v_d};
#last;
#last if($count == 250);
  }
}

sub blat_record {
  my ($self, $v_h) = @_;
  # now attempt the blat stuff
  my ($fh_target, $file_target) = tempfile( DIR => '/dev/shm', SUFFIX => '.fa', UNLINK => 1 );
  $self->blat_ref_alt($fh_target, $v_h);
  close $fh_target or die "Failed to close blat ref temp file";

  my $change_pos_low = $v_h->{change_pos};
  $change_pos_low++ if($v_h->{PC} eq 'I');
  my $range_l = ($v_h->{RE} - $v_h->{RS}) + 1;
  my $change_pos_high = $change_pos_low + $range_l; # REF based range, adjusted in func
  $v_h->{change_pos_low} = $change_pos_low;
  $v_h->{change_pos_high} = $change_pos_high;

  return $self->blat_reads($v_h, $file_target);
}

sub read_ranges {
  my ($self, $v_h) = @_;
  # return a string of chr:s-e... if approprate.
  my $ret_val;
  my ($chr, $q_start, $q_end) = ($v_h->{CHROM}, $v_h->{q_start}, $v_h->{q_end});
  my $read_buffer = $self->{max_insert};
  if($q_start + ($read_buffer * 2) > $q_end) {
    $ret_val = sprintf $LOCI_FMT, $chr, $q_start - $read_buffer, $q_end + $read_buffer;
  }
  else {
    $ret_val = sprintf $LOCI_FMT, $chr, $q_start - $read_buffer, $q_start + $read_buffer;
    $ret_val .= q{ };
    $ret_val .= sprintf $LOCI_FMT, $chr, $q_end - $read_buffer, $q_end + $read_buffer
  }
#print "$ret_val\n";
  return $ret_val;
}

sub blat_reads {
  my ($self, $v_h, $file_target) = @_;
  # setup ther temp files
  my ($fh_query, $file_query) = tempfile( DIR => '/dev/shm', SUFFIX => '.fa', UNLINK => 1);
  close $fh_query or die "Failed to close $file_query (query reads)";
  my ($fh_psl, $file_psl) = tempfile(DIR => '/dev/shm', SUFFIX => '.psl', UNLINK => 1);
  close $fh_psl or die "Failed to close $file_psl (psl output)";

  my $c_blat = sprintf $READS_AND_BLAT, $self->{hts_file}, $self->read_ranges($v_h), $file_query, $file_target, $file_query, $file_psl, $file_psl, $file_target, $file_query;
#print "$c_blat\n";
  my ($c_out, $c_err, $c_exit) = capture { system($c_blat); };
  if($c_exit) {
    warn "An error occurred while executing $c_blat\n";
    warn "\tERROR$c_err\n";
    exit $c_exit;
  }

# print "query.fa\n";
#system("cat $file_query");
# print "target.fa\n";
#system("cat $file_target");
#print "\n$c_out\n";
# exit 1;

  my ($wtp, $wtn, $mtp, $mtn) = $self->psl_axt_parser(\$c_out, $v_h);
  if($wtp + $wtn == 0) {
    ($wtp, $wtn) = $self->sam_depth($v_h);
  }
  my $mtr = $mtp+$mtn;
  my $depth = $wtp+$wtn+$mtr;
  my $vaf = $depth ? sprintf("%.3f", $mtr/$depth) : 0;
  return ($wtp, $wtn, $mtp, $mtn, $vaf);
}

sub psl_axt_parser {
  my ($self, $blat_axt, $v_h) = @_;
  # collate the data by readname and order by score
  my @lines = split /\n/, ${$blat_axt};
  my $line_c = @lines;
  # group all reads and order by score
  my %reads;
  for(my $i = 0; $i<$line_c; $i+=4) {
    my ($id, $t_name, $t_start, $t_end, $q_name, $q_start, $q_end, $strand, $score) = split q{ }, $lines[$i];
    push @{$reads{$q_name}{$score}}, [$t_name, $t_start, $t_end, $q_name, $q_start, $q_end, $strand, $score, $lines[$i+1], $lines[$i+2]];
  }

#print Dumper(\%reads);
my $SKIP_EVENT = q{};

  my %TYPE_STRAND;
  # sort keys for consistency
  READ: for my $read(sort keys %reads) {
    # get the alignment with the highest score
    for my $score(sort {$b<=>$a} keys %{$reads{$read}}) {
      my @records = @{$reads{$read}{$score}};
      if(@records != 1) { # if best score has more than one alignment it is irrelevant
#print Dumper(\@records);
        next READ;
      }
      my $record = $records[0];
      if($self->parse_axt_event($v_h, $record) == 1) {
        $TYPE_STRAND{$record->[0].$record->[6]} += 1;
      }
#$SKIP_EVENT = <STDIN> unless($SKIP_EVENT eq q{1});
#chomp $SKIP_EVENT;
      next READ; # remaining items are worse alignments
    }
  }
  my $wtp = $TYPE_STRAND{'REF+'} || 0;
  my $wtn = $TYPE_STRAND{'REF-'} || 0;
  my $mtp = $TYPE_STRAND{'ALT+'} || 0;
  my $mtn = $TYPE_STRAND{'ALT-'} || 0;
  return ($wtp, $wtn, $mtp, $mtn);
}

sub parse_axt_event {
  my ($self, $v_h, $rec) = @_;
  my ($t_name, $t_start, $t_end, $q_name, $q_start, $q_end, $strand, $score, $t_seq, $q_seq) = @{$rec};

  # specific to deletion class
  my $change_seq = $v_h->{change_ref};
  my $change_pos_high = $v_h->{change_pos_high};
  if($t_name eq 'ALT') {
    $change_pos_high -= $v_h->{LEN};
    $change_seq = $v_h->{change_alt};
  }

# eval {
#   my $exp_pos = $v_h->{change_pos_low} - $t_start; # duplicate of if block code
#   print join "\t", $t_name, $t_start, $t_end, $q_name, $q_start, $q_end, $strand, $score;
#   printf "\n%s\n%s\n", $t_seq, $q_seq;
#   print "exp_pos: $exp_pos\n";
#   printf "q_seqlen: %s\n", length $q_seq;
#   my $sub_q_seq = substr($q_seq, $exp_pos, length $change_seq);
#   my $lpad = q{ } x $exp_pos;
#   print $lpad.$change_seq."\n";
#   print $lpad.$sub_q_seq."\n";
#   my $boo = $t_start <= $v_h->{change_pos_low} && $t_end > $change_pos_high;
#   print "$boo: $t_start <= $v_h->{change_pos_low} && $t_end > $change_pos_high\n";

#   print "INDEX: ".index($q_seq, $change_seq, $exp_pos)."\n";
# }; print "OUT OF BOUNDS\n" if $@;


  # all the reads that don't span the range
  my $retval = 0;
  if($t_start <= $v_h->{change_pos_low} && $t_end > $change_pos_high) {
    # look for the change (or absence) where we expect it
    my $exp_pos = $v_h->{change_pos_low} - $t_start;
    if($exp_pos <= length $q_seq) {
      my $sub_q_seq = substr($q_seq, $exp_pos, length $change_seq);
      if(length $change_seq == length $sub_q_seq # same length
        && index($q_seq, q{-}) == -1 # no gaps
        && substr($change_seq,0,1) eq substr($sub_q_seq,0,1) # matching first base
        && substr($change_seq,-1,1) eq substr($sub_q_seq,-1,1) # matching last base
      ) {
        $retval = 1;
      }
    }
  }
#print "KEEP: $retval\n";
  return $retval;
}

sub sam_depth {
  my ($self, $v_h) = @_;
  my $mid_point = int ($v_h->{RS} + (($v_h->{RE} - $v_h->{RS})*0.5));
  my $read_search = sprintf $LOCI_FMT, $v_h->{CHROM}, $mid_point, $mid_point;
  my $c_samcount = sprintf $SAM_DEPTH_PN, $self->{hts_file}, $read_search;
  my ($c_out, $c_err, $c_exit) = capture { system($c_samcount); };
  if($c_exit) {
    warn "An error occurred while executing $c_samcount\n";
    warn "\tERROR$c_err\n";
    exit $c_exit;
  }
  return (split /\n/, $c_out);
}

sub flanking_ref {
  my ($self, $v_h) = @_;
  # use max_rl to extend before and after the range_start/end
  my $ref_left = $self->{fai}->get_sequence_no_length(
    sprintf $LOCI_FMT,
    $v_h->{CHROM},
    ($v_h->{POS} - $self->{max_rl})+1,
    $v_h->{POS},
  );
#print "$ref_left\n";

  my $ref_right = $self->{fai}->get_sequence_no_length(
    sprintf $LOCI_FMT,
    $v_h->{CHROM},
    $v_h->{END},
    $v_h->{END} + $self->{max_rl},
  );
#print "$ref_right\n";
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
  # THIS MAY NOT BE CORRECT NOW
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

# printf ">REF\n%s%s%s\n", $ref_left, $ref, $ref_right;
# printf ">ALT\n%s%s%s\n", $ref_left, $alt, $ref_right;
# printf "%s - %s - %s\n", $seq_left, $ref, $seq_right;
# printf "%s - %s - %s\n", $seq_left, $alt, $seq_right;
# printf "%d - %d = %d\n", $r_end, $r_start, $r_end - $r_start;
# printf "%s %s %s\n", $q_start, $q_end, $change_at;

  $v_h->{q_start} = $q_start;
  $v_h->{q_end} = $q_end;
  $v_h->{change_pos} = $change_at;
  $v_h->{change_ref} = $change_ref;
  $v_h->{change_alt} = $change_alt;

  return 1;
}
