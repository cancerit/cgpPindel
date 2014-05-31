package Sanger::CGP::PindelPostProcessing::FilterRules;

use strict;
use Tabix;
use Sanger::CGP::Pindel;

my %RULE_DESCS = ('F001' => { 'tag' =>'INFO/LEN',
                              'name' => 'F001',
                              'desc' => 'Pass if Mt > Wt Reads: Likely GERMLINE',
                              'test' => \&flag_001 },
                  'F002' => { 'tag'  => 'INFO/LEN',
                              'name' => 'F002',
                              'desc' => 'No Wt calls in variants over 4bp in length: Likely GERMLINE',
                              'test' => \&flag_002},
                  'F003' => { 'tag'  => 'INFO/LEN',
                              'name' => 'F003',
                              'desc' => 'Tum low call count strand bias check',
                              'test' => \&flag_003},
                  'F004' => { 'tag' =>'INFO/LEN',
                              'name' => 'F004',
                              'desc' => 'Tum medium read depth strand bias check: Calls In 8% Reads Bt Depth 10 And 200 (inclusive)',
                              'test' => \&flag_004 },
                  'F005' => { 'tag'  => 'INFO/LEN',
                              'name' => 'F005',
                              'desc' => 'Tum high read depth strand bias check: Calls In 4% Reads > Depth 200',
                              'test' => \&flag_005},
                  'F006' => { 'tag'  => 'INFO/LEN',
                              'name' => 'F006',
                              'desc' => 'Small call excessive repeat check: Fail if Length <= 4 and Repeats > 9',
                              'test' => \&flag_006},
                  'F007' => { 'tag'  => 'INFO/LEN',
                              'name' => 'F007',
                              'desc' => 'Sufficient Normal Depth: If Mt Depth > 5 then Wt > 8% tum depth',
                              'test' => \&flag_007},
                  'F008' => { 'tag'  => 'INFO/REP',
                              'name' => 'F008',
                              'desc' => 'Wildtype contamination: Fail when wt reads > 5% mt reads.',
                              'test' => \&flag_008},
                  'F009' => { 'tag'  => 'INFO/LEN',
                              'name' => 'F009',
                              'desc' => 'Is coding: Pass when in gene footprint.',
                              'test' => \&flag_009},
                  'F010' => { 'tag'  => 'INFO/LEN',
                              'name' => 'F010',
                              'desc' => 'Variant must not exist within the Unmatched Normal Panel',
                              'test' => \&flag_010},
                  # 11 legacy
                  'F012' => { 'tag'  => 'INFO/LEN',
                              'name' => 'F012',
                              'desc' => 'Germline: When length < 11 and depth > 9, fail if the variant is seen in both 20% of normal reads AND 20% of tumour reads in either pindel or bwa',
                              'test' => \&flag_012},
                  # 13,14 legacy
                  'F015' => { 'tag'  => 'INFO/LEN',
                              'name' => 'F015',
                              'desc' => 'No normal calls',
                              'test' => \&flag_015},
                  'F016' => { 'tag'  => 'INFO/REP',
                              'name' => 'F016',
                              'desc' => 'Verify indel condition',
                              'test' => \&flag_016},
                  'F017' => { 'tag'  => 'INFO/LEN',
                              'name' => 'F017',
                              'desc' => 'Variant must not overlap with a simple repeat',
                              'test' => \&flag_017},
                  'F018' => { 'tag'  => 'INFO/LEN',
                              'name' => 'F018',
                              'desc' => 'Sufficient Depth: Pass if depth > 10',
                              'test' => \&flag_018},
);

our $previous_format_hash;
our $previous_format_string = q{};
our $vcf_flagging_repeats_tabix;
our $vcf_flagging_unmatched_normals_tabix;

sub rule {
  my $rule = shift;
  return $RULE_DESCS{$rule};
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
    $vcf_flagging_unmatched_normals_tabix = new Tabix(-data => $ENV{VCF_FLAGGING_UNMATCHED_NORMALS},-index => $ENV{VCF_FLAGGING_UNMATCHED_NORMALS}.'.tbi');
  }
}

sub reuse_repeats_tabix {
  unless(defined $vcf_flagging_repeats_tabix) {
    $vcf_flagging_repeats_tabix = new Tabix(-data => $ENV{VCF_FLAGGING_REPEATS},-index => $ENV{VCF_FLAGGING_REPEATS}.'.tbi');
  }
}

sub flag_001 {
  my ($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF) = @_;
  use_prev($$RECORD[8]);

  my @nor_geno = split(':',$$RECORD[9]);
  my @tum_geno = split(':',$$RECORD[10]);

  if($nor_geno[$previous_format_hash->{'PP'}] + $nor_geno[$previous_format_hash->{'NP'}]  >=
    $tum_geno[$previous_format_hash->{'PP'}] + $tum_geno[$previous_format_hash->{'NP'}]){
    return $FAIL;
  }
  return $PASS;
}

sub flag_002 {
  my ($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF) = @_;
  return $PASS if($MATCH <= 4);
  use_prev($$RECORD[8]);

  my @nor_geno = split(':',$$RECORD[9]);

  if($nor_geno[$previous_format_hash->{'PP'}] + $nor_geno[$previous_format_hash->{'NP'}] > 0) {
    return $FAIL;
  }
  return $PASS;
}

sub flag_003 {
  my ($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF) = @_;
  use_prev($$RECORD[8]);

  my @tum_geno = split(':',$$RECORD[10]);
  if(($tum_geno[$previous_format_hash->{'PP'}] >= 3 || $tum_geno[$previous_format_hash->{'NP'}] >= 3)){
    return $PASS;
  }

  if(($tum_geno[$previous_format_hash->{'PP'}] >= 2 && $tum_geno[$previous_format_hash->{'NP'}] >= 2)){
    return $PASS;
  }

  return $FAIL;
}

sub flag_004 {
  my ($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF) = @_;
  use_prev($$RECORD[8]);

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

  my $nor_d = $nor_geno[$previous_format_hash->{'PR'}] + $nor_geno[$previous_format_hash->{'NR'}];
  my $tum_d = $tum_geno[$previous_format_hash->{'PR'}] + $tum_geno[$previous_format_hash->{'NR'}];

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

  if($nor_geno[$previous_format_hash->{'PP'}] + $nor_geno[$previous_format_hash->{'NP'}] >
    ($tum_geno[$previous_format_hash->{'PP'}] + $tum_geno[$previous_format_hash->{'NP'}]) * 0.05){
      return $FAIL;

  }
  return $PASS;
}

sub flag_009 {
  my ($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF) = @_;
  ### HACK Dirty dirty dirty......
  unless($main::VCF_IS_CODING_TABIX){
    use Tabix;
    $main::VCF_IS_CODING_TABIX = new Tabix(-data => $ENV{VCF_IS_CODING},-index => $ENV{VCF_IS_CODING}.'.tbi');
  }

  my $ret = eval{
    ## half open interval based.... i.e. from is zero based to is one based.
    # as vcf POS for indels is the previous base pos is already 0-based
    my $res = $main::VCF_IS_CODING_TABIX->query($CHROM,$POS,($POS+$MATCH));
    return $FAIL if(!defined $res->get); # no valid entries (chromosome not in index) so must FAIL
    return $PASS if($main::VCF_IS_CODING_TABIX->read($res));
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
  $from -= ($length_off + 1); # additional -1 to switch to 0 based
  $to += $length_off;

  my $ret = eval{
    my $res = $vcf_flagging_unmatched_normals_tabix->query($CHROM,$from,$to);

    return $PASS if(!defined $res->get); # no valid entries (chromosome not in index) so must pass

    while(my $line = $vcf_flagging_unmatched_normals_tabix->read($res)){
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

      my $tum_target_depth = ($tum_geno[$previous_format_hash->{'ND'}] + $tum_geno[$previous_format_hash->{'PD'}]) * 0.2;
      my $nor_target_depth = ($nor_geno[$previous_format_hash->{'ND'}] + $nor_geno[$previous_format_hash->{'PD'}]) * 0.2;

      if((($nor_geno[$previous_format_hash->{'NP'}] + $nor_geno[$previous_format_hash->{'PP'}] >= $nor_target_depth) &&
        ($tum_geno[$previous_format_hash->{'NP'}] + $tum_geno[$previous_format_hash->{'PP'}] >= $tum_target_depth))
        ||
        (($nor_geno[$previous_format_hash->{'NB'}] + $nor_geno[$previous_format_hash->{'PB'}] >= $nor_target_depth) &&
        ($tum_geno[$previous_format_hash->{'NB'}] + $tum_geno[$previous_format_hash->{'PB'}] >= $tum_target_depth))
        ){
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

  if($nor_geno[$previous_format_hash->{'PU'}] + $nor_geno[$previous_format_hash->{'NU'}] == 0){
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
    my $res = $vcf_flagging_repeats_tabix->query($CHROM,($from-1),$to);
    return $PASS if(!defined $res->get); # no valid entries (chromosome not in index) so must pass
    return $FAIL if($vcf_flagging_repeats_tabix->read($res));
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

1;
