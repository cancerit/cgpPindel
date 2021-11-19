package Sanger::CGP::PindelPostProcessing::FragmentFilterRules;

########## LICENCE ##########
# Copyright (c) 2014-2021 Genome Research Ltd.
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
use Bio::DB::HTS::Tabix;
use Sanger::CGP::Pindel;

my %RULE_DESCS = ('FF001' => { 'tag' =>'INFO/LEN',
                              'name' => 'FF001',
                              'desc' => 'Pass if Mt > Wt Reads: Likely GERMLINE',
                              'test' => \&flag_001 },
                  'FF002' => { 'tag'  => 'INFO/LEN',
                              'name' => 'FF002',
                              'desc' => 'No Wt calls in variants over 4bp in length: Likely GERMLINE',
                              'test' => \&flag_002},
                  'FF003' => { 'tag'  => 'INFO/LEN',
                              'name' => 'FF003',
                              'desc' => 'Tum low call count strand bias check',
                              'test' => \&flag_003},
                  'FF004' => { 'tag' =>'INFO/LEN',
                              'name' => 'FF004',
                              'desc' => 'Tum medium read depth strand bias check: Calls In 8% Reads Bt Depth 10 And 200 (inclusive)',
                              'test' => \&flag_004 },
                  'FF005' => { 'tag'  => 'INFO/LEN',
                              'name' => 'FF005',
                              'desc' => 'Tum high read depth strand bias check: Calls In 4% Reads > Depth 200',
                              'test' => \&flag_005},
                  'FF006' => { 'tag'  => 'INFO/LEN',
                              'name' => 'FF006',
                              'desc' => 'Small call excessive repeat check: Fail if Length <= 4 and Repeats > 9',
                              'test' => \&flag_006},
                  'FF007' => { 'tag'  => 'INFO/LEN',
                              'name' => 'FF007',
                              'desc' => 'Sufficient Normal Depth: If Mt Depth > 5 then Wt > 8% tum depth',
                              'test' => \&flag_007},
                  'FF008' => { 'tag'  => 'INFO/REP',
                              'name' => 'FF008',
                              'desc' => 'Wildtype contamination: Fail when wt reads > 5% mt reads.',
                              'test' => \&flag_008},
                  'FF009' => { 'tag'  => 'INFO/LEN',
                              'name' => 'FF009',
                              'desc' => 'Is coding: Pass when in gene footprint.',
                              'test' => \&flag_009},
                  'FF010' => { 'tag'  => 'INFO/LEN',
                              'name' => 'FF010',
                              'desc' => 'Variant must not exist within the Unmatched Normal Panel',
                              'test' => \&flag_010},
                  # 11 legacy
                  'FF012' => { 'tag'  => 'INFO/LEN',
                              'name' => 'FF012',
                              'desc' => 'Germline: When length < 11 and depth > 9, fail if the variant is seen in both 20% of normal reads AND 20% of tumour reads in either pindel or bwa',
                              'test' => \&flag_012},
                  # 13,14 legacy
                  'FF015' => { 'tag'  => 'INFO/LEN',
                              'name' => 'FF015',
                              'desc' => 'No normal calls',
                              'test' => \&flag_015},
                  'FF016' => { 'tag'  => 'INFO/REP',
                              'name' => 'FF016',
                              'desc' => 'Verify indel condition',
                              'test' => \&flag_016},
                  'FF017' => { 'tag'  => 'INFO/LEN',
                              'name' => 'FF017',
                              'desc' => 'Variant must not overlap with a simple repeat',
                              'test' => \&flag_017},
                  'FF018' => { 'tag'  => 'INFO/LEN',
                              'name' => 'FF018',
                              'desc' => 'Sufficient Depth: Pass if depth > 10',
                              'test' => \&flag_018},
                  'FF019' => { 'tag'  => 'INFO/LEN',
                              'name' => 'FF019',
                              'desc' => 'Fail when tumour supporting fragments < 3 or tumour fraction of supporting fragments < 0.05',
                              'test' => \&flag_019},
                    'FF020' => { 'tag'  => 'INFO/LEN',
                              'name' => 'FF020',
                              'desc' => 'Allow some contamination in matched normal due to FFPR block acquired samples and allow for low level sequencing/PCR artefacts',
                              'test' => \&flag_020},
);

our $previous_format_hash;
our $previous_format_string = q{};
our $vcf_flagging_repeats_tabix;
our $vcf_flagging_unmatched_normals_tabix;

sub rule {
  my (undef, $rule) = @_; # being called like an object function so throw away first varaible
  return $RULE_DESCS{$rule};
}

sub available_rules {
  return sort keys %RULE_DESCS;
}

sub use_prev {
  my $format = shift;
  ### HACK Dirty dirty dirty...... done to try and cut down the number of times I have to parse the FORMAT string I am storing it as a global variable.
  if($format ne $previous_format_string){
    my $i = 0;
    map {$previous_format_hash->{$_} = $i++} split(':',$format);
    $previous_format_string = $format;
  }
}

sub reuse_unmatched_normals_tabix {
  unless(defined $vcf_flagging_unmatched_normals_tabix){
    $vcf_flagging_unmatched_normals_tabix = new Bio::DB::HTS::Tabix(filename=> $ENV{VCF_FLAGGING_UNMATCHED_NORMALS});
  }
}

sub reuse_repeats_tabix {
  unless(defined $vcf_flagging_repeats_tabix) {
    $vcf_flagging_repeats_tabix = new Bio::DB::HTS::Tabix(filename=> $ENV{VCF_FLAGGING_REPEATS});
  }
}

sub flag_001 {
  my ($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF) = @_;
  use_prev($$RECORD[8]);

  my @nor_geno = split(':',$$RECORD[9]);
  my @tum_geno = split(':',$$RECORD[10]);

  if($nor_geno[$previous_format_hash->{'FC'}] >= $tum_geno[$previous_format_hash->{'FC'}]){
    return $FAIL;
  }
  return $PASS;
}

sub flag_002 {
  my ($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF) = @_;
  return $PASS if($MATCH <= 4);
  use_prev($$RECORD[8]);

  my @nor_geno = split(':',$$RECORD[9]);

  if($nor_geno[$previous_format_hash->{'FC'}] > 0) {
    return $FAIL;
  }
  return $PASS;
}

sub flag_003 {
  my ($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF) = @_;
  use_prev($$RECORD[8]);

  my @tum_geno = split(':',$$RECORD[10]);
  if(($tum_geno[$previous_format_hash->{'PU'}] >= 3 || $tum_geno[$previous_format_hash->{'NU'}] >= 3)){
    return $PASS;
  }

  if(($tum_geno[$previous_format_hash->{'PU'}] >= 2 && $tum_geno[$previous_format_hash->{'NU'}] >= 2)){
    return $PASS;
  }

  return $FAIL;
}

sub flag_004 {
  my ($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF) = @_;
  use_prev($$RECORD[8]);

  # unchanged, not fragment based
  my $ret = $FAIL;
  my @tum_geno = split(':',$$RECORD[10]);

  my $prd = $tum_geno[$previous_format_hash->{'PR'}];
  my $nrd = $tum_geno[$previous_format_hash->{'NR'}];
  my $rd = $prd + $nrd;
  my $pp = $tum_geno[$previous_format_hash->{'PU'}];
  my $np = $tum_geno[$previous_format_hash->{'NU'}];
  my $p = $pp + $np;

  if($rd >= 10 && $rd < 200) {
    if($pp > 0 && $np > 0) {
      $ret = $PASS if($p >= ($rd * 0.05));
    }
    elsif($pp > 0) {
      if($prd) {
        $ret = $PASS if($pp >= ($prd * 0.08));
      }
      else {
        $ret = $PASS;
      }
    }
    elsif($np > 0){
      if($nrd) {
        $ret = $PASS if($np >= ($nrd * 0.08));
      }
      else {
        $ret = $PASS;
      }
    }
  }
  else {
    $ret = $PASS;
  }
  return $ret;
}

sub flag_005 {
  my ($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF) = @_;
  use_prev($$RECORD[8]);

  # unchanged, not fragment based
  my $ret = $FAIL;
  my @tum_geno = split(':',$$RECORD[10]);

  my $prd = $tum_geno[$previous_format_hash->{'PR'}];
  my $nrd = $tum_geno[$previous_format_hash->{'NR'}];
  my $rd = $prd + $nrd;
  my $pp = $tum_geno[$previous_format_hash->{'PU'}];
  my $np = $tum_geno[$previous_format_hash->{'NU'}];
  my $p = $pp + $np;

  if($rd >= 200){
    if($pp > 0 && $np > 0){
      $ret = $PASS if($p >= ($rd * 0.04));
    }
    elsif($pp > 0){
      if($prd){
        $ret = $PASS if($pp >= ($prd * 0.04));
      }
      else{
        $ret = $PASS;
      }
    }
    elsif($np > 0){
      if($nrd){
        $ret = $PASS if($np >= ($nrd * 0.04));
      }
      else{
        $ret = $PASS;
      }
    }
  }
  else{
    $ret = $PASS;
  }
  return $ret;
}

sub flag_006 {
  my ($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF) = @_;
  if($MATCH <= 4){
    my ($rep) = $$RECORD[7] =~ /REP=(\d+)/;
    if($rep > 9) {
      return $FAIL;
    }
  }
  return $PASS;
}

sub flag_007 {
  my ($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF) = @_;
  use_prev($$RECORD[8]);

  my @nor_geno = split(':',$$RECORD[9]);
  my @tum_geno = split(':',$$RECORD[10]);

  my $nor_d = $nor_geno[$previous_format_hash->{'FD'}];
  my $tum_d = $tum_geno[$previous_format_hash->{'FD'}];

  if($tum_d > 5){
    if($nor_d >= ($tum_d * 0.08)){
      return $PASS;
    }
  }

  return $FAIL;
}

sub flag_008 {
  my ($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF) = @_;
  use_prev($$RECORD[8]);

  my @nor_geno = split(':',$$RECORD[9]);
  my @tum_geno = split(':',$$RECORD[10]);

  if($nor_geno[$previous_format_hash->{'FC'}] > ($tum_geno[$previous_format_hash->{'FC'}] * 0.05)){
      return $FAIL;

  }
  return $PASS;
}

sub flag_009 {
  my ($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF) = @_;
  ### HACK Dirty dirty dirty......
  unless($main::VCF_IS_CODING_TABIX){
    use Bio::DB::HTS::Tabix;
    $main::VCF_IS_CODING_TABIX = new Bio::DB::HTS::Tabix(filename=> $ENV{VCF_IS_CODING});
  }

  my $ret = eval{
    # as vcf POS for indels is the previous base pos is 0-based, but the new TABIX requires 1-based
    my $iter = $main::VCF_IS_CODING_TABIX->query_full($CHROM,$POS+1,($POS+$MATCH));
    return $FAIL if(!defined $iter); # no valid entries (chromosome not in index) so must FAIL
    while($iter->next){
      return $PASS;
    }
    return $FAIL;
  };
  if($@) {
    die $@;
  }
  return $ret;
}

sub flag_010 {
  my ($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF) = @_;
  reuse_unmatched_normals_tabix();

  my $length_off = ($MATCH <= 2) ? 1 : 20;

  my ($from) = ";$$RECORD[7]" =~ m/;RS=([^;]+)/;
  my ($to) = ";$$RECORD[7]" =~ m/;RE=([^;]+)/;
  #Range is the bases surrounding a position so need to bump back to the actual repetitive tract (or single pos)
  # but only when REP is > 0
  if(";$$RECORD[7];" !~ m/;REP=0;/) {
    $from++;
    $to--;
  }
  # then apply the range fudging
  $from -= $length_off; # no fudge of position as tabix for bed is now 1 based
  $to += $length_off;

  my $ret = eval{
    my $iter = $vcf_flagging_unmatched_normals_tabix->query_full($CHROM,$from,$to);
    return $PASS if(!defined $iter); # no valid entries (chromosome not in index) so must pass
    while($iter->next){
      return $FAIL;
    }
    return $PASS;
  };
  if($@) {
    die $@;
  }
  return $ret;
}

sub flag_012 {
  my ($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF) = @_;
  if($MATCH < 11){
    use_prev($$RECORD[8]);

    my @nor_geno = split(':',$$RECORD[9]);
    my @tum_geno = split(':',$$RECORD[10]);

    if(($tum_geno[$previous_format_hash->{'ND'}] + $tum_geno[$previous_format_hash->{'PD'}]) > 9 &&
       ($nor_geno[$previous_format_hash->{'ND'}] + $nor_geno[$previous_format_hash->{'PD'}]) > 9){

      my $tum_target_depth = $tum_geno[$previous_format_hash->{'FD'}] * 0.2;
      my $nor_target_depth = $nor_geno[$previous_format_hash->{'FD'}]* 0.2;

      if(($nor_geno[$previous_format_hash->{'FC'}] >= $nor_target_depth) &&
         ($tum_geno[$previous_format_hash->{'FC'}] >= $tum_target_depth)){
        return $FAIL;
      }
    }

  }
  return $PASS;
}

sub flag_015 {
  my ($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF) = @_;
  use_prev($$RECORD[8]);

  my @nor_geno = split(':',$$RECORD[9]);

  if($nor_geno[$previous_format_hash->{'FC'}] == 0){
    return $PASS;
  }
  return $FAIL;
}

sub flag_016 {
  my ($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF) = @_;
  use_prev($$RECORD[8]);

  my @nor_geno = split(':',$$RECORD[9]);
  my @tum_geno = split(':',$$RECORD[10]);

  if($tum_geno[$previous_format_hash->{'PP'}] + $tum_geno[$previous_format_hash->{'NP'}] > 4){
    if($tum_geno[$previous_format_hash->{'PB'}] + $tum_geno[$previous_format_hash->{'NB'}] > 0){
      return $PASS;
    }else{

      ## RECORD[3]: ref
      ## RECORD[4]: alt
      my $comp = length($$RECORD[3]) > 1 && length($$RECORD[4]) == 1 ? 1 : 0;

      if($MATCH == $comp){ ##repeats...

        if($tum_geno[$previous_format_hash->{'PP'}] > 0 &&
        $tum_geno[$previous_format_hash->{'NP'}] > 0){
          return $PASS;
        }
      }
    }
  }
  return $FAIL;
}

sub flag_017 {
  my ($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF) = @_;
  reuse_repeats_tabix();

  my($from) = ";$$RECORD[7]" =~ m/;RS=(\d+)/;
  my($to) = ";$$RECORD[7]" =~ m/;RE=(\d+)/;

  my $ret = eval{
    # as vcf POS for indels is the previous base pos is 0-based, but the new TABIX requires 1-based
    my $iter = $vcf_flagging_repeats_tabix->query_full($CHROM,$from,$to);
    return $PASS if(!defined $iter); # no valid entries (chromosome not in index) so must pass
    while($iter->next){
      return $FAIL;
    }
    return $PASS;
  };
  if($@) {
    die $@;
  }
  return $ret;
}

sub flag_018 {
  my ($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF) = @_;
  use_prev($$RECORD[8]);

  my @nor_geno = split(':',$$RECORD[9]);
  my @tum_geno = split(':',$$RECORD[10]);
  if(($nor_geno[$previous_format_hash->{'PR'}] + $nor_geno[$previous_format_hash->{'NR'}] >= 10) &&
    ($tum_geno[$previous_format_hash->{'PR'}] + $tum_geno[$previous_format_hash->{'NR'}] >= 10)){
    return $PASS;
  }

  return $FAIL;
}

sub flag_019 {
  my ($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF) = @_;
  use_prev($$RECORD[8]);

  my @tum_geno = split(':',$$RECORD[10]);

  if($tum_geno[$previous_format_hash->{'FC'}] < 3){
    return $FAIL;
  }
  # previous test confirms FC/FD can't be 0, so no div0 check required
  if ($tum_geno[$previous_format_hash->{'FC'}] / $tum_geno[$previous_format_hash->{'FD'}] < 0.05){
    return $FAIL;
  }

  return $PASS;
}

sub flag_020 {
  my ($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF) = @_;
  use_prev($$RECORD[8]);

  my @nor_geno = split(':',$$RECORD[9]);
  my @tum_geno = split(':',$$RECORD[10]);

  my $fd_total = $nor_geno[$previous_format_hash->{'FD'}] + $tum_geno[$previous_format_hash->{'FD'}];

  if($fd_total < 200 &&
    $nor_geno[$previous_format_hash->{'FC'}] <= 1 && 
    $nor_geno[$previous_format_hash->{'FD'}] >= 10 && 
    $nor_geno[$previous_format_hash->{'FC'}] < ($tum_geno[$previous_format_hash->{'FC'}] * 0.1)
  ){
    return $FAIL;
  }
  
  my $tumfc_over_tumfd = $tum_geno[$previous_format_hash->{'FD'}] > 0 ? $tum_geno[$previous_format_hash->{'FC'}] / $tum_geno[$previous_format_hash->{'FD'}] : undef;
  my $norfc_over_norfd = $nor_geno[$previous_format_hash->{'FD'}] > 0 ? $nor_geno[$previous_format_hash->{'FC'}] / $nor_geno[$previous_format_hash->{'FD'}] : undef;

  if($fd_total < 200){
    if(($nor_geno[$previous_format_hash->{'FC'}] == 1 || $nor_geno[$previous_format_hash->{'FC'}] == 2) &&
      $norfc_over_norfd <= 0.05 &&
      $tumfc_over_tumfd >= 0.2
    ){
    return $FAIL;
    }
    
  }else{
    if($norfc_over_norfd > 0.02 && $tumfc_over_tumfd < 0.2){
      return $FAIL;
    }
  }

  return $PASS;
}

1;
