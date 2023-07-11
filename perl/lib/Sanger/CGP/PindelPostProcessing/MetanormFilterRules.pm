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
package Sanger::CGP::PindelPostProcessing::MetanormFilterRules;

use strict;
use Bio::DB::HTS::Tabix;
use Sanger::CGP::Pindel;

my %RULE_DESCS = (
                  'F006' => { 'tag'  => 'INFO/LEN',
                              'name' => 'F006',
                              'desc' => 'Small call excessive repeat check: Fail if Length <= 4 and Repeats > 9',
                              'test' => \&flag_006},
                  'F017' => { 'tag'  => 'INFO/LEN',
                              'name' => 'F017',
                              'desc' => 'Variant must not overlap with a simple repeat',
                              'test' => \&flag_017},
                  'LONG' => { 'tag'  => 'INFO/LEN',
                              'name' => 'LONG',
                              'desc' => 'Event larger than 1kbp',
                              'test' => \&long}

);

our $previous_format_hash;
our $previous_format_string = q{};
our $vcf_flagging_repeats_tabix;

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

sub reuse_repeats_tabix {
  unless(defined $vcf_flagging_repeats_tabix) {
    $vcf_flagging_repeats_tabix = new Bio::DB::HTS::Tabix(filename=> $ENV{VCF_FLAGGING_REPEATS});
  }
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

sub long {
  my ($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF) = @_;
  if($MATCH > 1000 ){
    return $FAIL
  }
  return $PASS;
}


1;
