####################################################
# Copyright (c) 2014-2021 Genome Research Ltd.
# Author: CASM/Cancer IT, cgphelp@sanger.ac.uk
# See LICENCE for details
####################################################

use strict;
use warnings;
use Cwd 'abs_path';
use Test::More;
use Data::Dumper;
use Const::Fast qw(const);

const my @AVAILABLE_RULES => qw(FF001 FF002 FF003 FF004 FF005 FF006 FF007 FF008 FF009 FF010 FF012 FF015 FF016 FF017 FF018 FF019 FF020);

my %rule_test_dispatch = ('FF001' => \&_test_FF001,
                          'FF002' => \&_test_FF002,
                          'FF003' => \&_test_FF003,
                          'FF004' => \&_test_FF004,
                          'FF005' => \&_test_FF005,
                          'FF006' => \&_test_FF006,
                          'FF007' => \&_test_FF007,
                          'FF008' => \&_test_FF008,
                          'FF009' => \&_test_FF009,
                          'FF010' => \&_test_FF010,
                          'FF012' => \&_test_FF012,
                          'FF015' => \&_test_FF015,
                          'FF016' => \&_test_FF016,
                          'FF017' => \&_test_FF017,
                          'FF018' => \&_test_FF018,
                          'FF019' => \&_test_FF019,
                          'FF020' => \&_test_FF020,
                        );

use_ok('Sanger::CGP::PindelPostProcessing::VcfSoftFlagger');
use_ok('Sanger::CGP::PindelPostProcessing::FragmentFilterRules');

my @rules_found = Sanger::CGP::PindelPostProcessing::FragmentFilterRules::available_rules();
is_deeply(\@rules_found, \@AVAILABLE_RULES, 'Expected set of rules are implemented');
for my $flag(@AVAILABLE_RULES) {
  $rule_test_dispatch{$flag}(Sanger::CGP::PindelPostProcessing::FragmentFilterRules->rule($flag));
}

done_testing();

sub _test_FF001{
	my ($filter_hash) = @_;
  subtest "Test rule FF001" => sub {
    my $test1 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	FC	5	5';
    my $test2 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	FC	6	4';
    my $test3 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	FC	4	6';
    my $sub = $filter_hash->{test};

    my($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF);
    $CHROM  = [split("\t",$test1)]->[0];
    $POS    = [split("\t",$test1)]->[1];
    $FAIL   = 1;
    $PASS   = 0;
    #$VCF    = $$opts{vcf};

    is($filter_hash->{tag}, 'INFO/LEN',"_test_FF001 check the correct info tag has been set for the rule");

    $RECORD = [split("\t",$test1)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_FF001 pWt == pMt");
    $RECORD = [split("\t",$test2)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_FF001 pWt > pMt");
    $RECORD = [split("\t",$test3)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_FF001 pWt < pMt");
	};
}


sub _test_FF002{
	my ($filter_hash) = @_;
  subtest "Test rule FF002" => sub {
    my $test1 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=4;SM=138;S1=10;S2=203.791;REP=18	FC	1	0';
    my $test2 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=4;SM=138;S1=10;S2=203.791;REP=18	FC	0	0';

    my $test3 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=5;SM=138;S1=10;S2=203.791;REP=18	FC	1	0';
    my $test4 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=5;SM=138;S1=10;S2=203.791;REP=18	FC	0	0';


    my $sub = $filter_hash->{test};

    my($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF);
    $CHROM  = [split("\t",$test1)]->[0];
    $POS    = [split("\t",$test1)]->[1];
    $FAIL   = 1;
    $PASS   = 0;
    #$VCF    = $$opts{vcf};

    is($filter_hash->{tag}, 'INFO/LEN',"_test_FF002 check the correct info tag has been set for the rule");

    $RECORD = [split("\t",$test1)];
    $MATCH = 4;
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_FF002 length == 4 pWt == 1");
    $RECORD = [split("\t",$test2)];
    $MATCH = 4;
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_FF002 length == 4 pWt == 0");

    $RECORD = [split("\t",$test3)];
    $MATCH = 5;
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_FF002 length == 5 pWt == 1");
    $RECORD = [split("\t",$test4)];
    $MATCH = 5;
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_FF002 length == 5 pWt == 0");
  }
}

sub _test_FF003{
	my ($filter_hash) = @_;
  subtest "Test rule FF003" => sub {
    my $test1 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PU:NU	0:0	3:0';
    my $test2 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PU:NU	0:0	0:3';
    my $test3 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PU:NU	0:0	4:0';
    my $test4 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PU:NU	0:0	0:4';
    my $test5 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PU:NU	0:0	3:4';
    my $test6 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PU:NU	0:0	2:2';
    my $test7 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PU:NU	0:0	2:1';
    my $test8 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PU:NU	0:0	1:2';
    my $test9 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PU:NU	0:0	2:0';
    my $test10 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PU:NU	0:0	0:2';
    my $test11 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PU:NU	0:0	0:0';



    my $sub = $filter_hash->{test};

    my($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF);
    $CHROM  = [split("\t",$test1)]->[0];
    $POS    = [split("\t",$test1)]->[1];
    $FAIL   = 1;
    $PASS   = 0;
    #$VCF    = $$opts{vcf};

    is($filter_hash->{tag}, 'INFO/LEN',"_test_FF003 check the correct info tag has been set for the rule");

    $RECORD = [split("\t",$test1)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_FF003 pMtPos == 3 pMtNeg == 0");
    $RECORD = [split("\t",$test2)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_FF003 pMtPos == 0 pMtNeg == 3");
    $RECORD = [split("\t",$test3)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_FF003 pMtPos == 4 pMtNeg == 0");
    $RECORD = [split("\t",$test4)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_FF003 pMtPos == 0 pMtNeg == 4");

    $RECORD = [split("\t",$test5)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_FF003 pMtPos == 3 pMtNeg == 4");
    $RECORD = [split("\t",$test6)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_FF003 pMtPos == 2 pMtNeg == 2");
    $RECORD = [split("\t",$test7)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_FF003 pMtPos == 2 pMtNeg == 1");
    $RECORD = [split("\t",$test8)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_FF003 pMtPos == 1 pMtNeg == 2");
    $RECORD = [split("\t",$test9)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_FF003 pMtPos == 2 pMtNeg == 0");
    $RECORD = [split("\t",$test10)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_FF003 pMtPos == 0 pMtNeg == 2");
    $RECORD = [split("\t",$test11)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_FF003 pMtPos == 0 pMtNeg == 0");
	};
}


sub _test_FF004{
	my ($filter_hash) = @_;
  subtest "Test rule FF004" => sub {
    my $test1 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PR:NR:PU:NU	0:0:0:0	9:0:0:0';
    my $test2 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PR:NR:PU:NU	0:0:0:0	10:0:0:0';
    my $test3 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PR:NR:PU:NU	0:0:0:0	199:0:0:0';
    my $test4 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PR:NR:PU:NU	0:0:0:0	200:0:0:0';

    my $test5 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PR:NR:PU:NU	0:0:0:0	100:0:2:3';
    my $test6 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PR:NR:PU:NU	0:0:0:0	100:0:2:2';


    my $test7 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PR:NR:PU:NU	0:0:0:0	0:100:4:0';
    my $test8 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PR:NR:PU:NU	0:0:0:0	100:0:8:0';
    my $test9 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PR:NR:PU:NU	0:0:0:0	100:0:7:0';


    my $test10 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PR:NR:PU:NU	0:0:0:0	100:0:0:4';
    my $test11 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PR:NR:PU:NU	0:0:0:0	0:100:0:8';
    my $test12 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PR:NR:PU:NU	0:0:0:0	0:100:0:7';


    my $sub = $filter_hash->{test};

    my($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF);
    $CHROM  = [split("\t",$test1)]->[0];
    $POS    = [split("\t",$test1)]->[1];
    $FAIL   = 1;
    $PASS   = 0;
    #$VCF    = $$opts{vcf};

    is($filter_hash->{tag}, 'INFO/LEN',"_test_FF004 check the correct info tag has been set for the rule");

    $RECORD = [split("\t",$test1)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_FF004 ");
    $RECORD = [split("\t",$test2)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_FF004 ");
    $RECORD = [split("\t",$test3)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_FF004 ");
    $RECORD = [split("\t",$test4)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_FF004 ");

    $RECORD = [split("\t",$test5)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_FF004 ");
    $RECORD = [split("\t",$test6)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_FF004 ");

    $RECORD = [split("\t",$test7)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_FF004 ");

    $RECORD = [split("\t",$test8)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_FF004 ");
    $RECORD = [split("\t",$test9)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_FF004 ");

    $RECORD = [split("\t",$test10)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_FF004 ");

    $RECORD = [split("\t",$test11)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_FF004 ");
    $RECORD = [split("\t",$test12)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_FF004 ");
  };
}

sub _test_FF005{
	my ($filter_hash) = @_;
  subtest "Test rule FF005" => sub {
    my $test1 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PR:NR:PU:NU	0:0:0:0	199:0:1:0';
    my $test2 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PR:NR:PU:NU	0:0:0:0	0:199:1:0';

    my $test3 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PR:NR:PU:NU	0:0:0:0	200:0:4:4';
    my $test4 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PR:NR:PU:NU	0:0:0:0	0:200:4:4';
    my $test5 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PR:NR:PU:NU	0:0:0:0	199:1:4:4';
    my $test6 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PR:NR:PU:NU	0:0:0:0	1:199:4:4';

    my $test7 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PR:NR:PU:NU	0:0:0:0	200:0:1:1';
    my $test8 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PR:NR:PU:NU	0:0:0:0	200:0:3:4';
    my $test9 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PR:NR:PU:NU	0:0:0:0	200:0:4:4';
    my $test10 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PR:NR:PU:NU	0:0:0:0	200:0:4:5';

    my $test11 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PR:NR:PU:NU	0:0:0:0	0:200:1:0';
    my $test12 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PR:NR:PU:NU	0:0:0:0	200:0:7:0';
    my $test13 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PR:NR:PU:NU	0:0:0:0	200:0:8:0';

    my $test14 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PR:NR:PU:NU	0:0:0:0	200:0:0:1';
    my $test15 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PR:NR:PU:NU	0:0:0:0	0:200:0:7';
    my $test16 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PR:NR:PU:NU	0:0:0:0	0:200:0:8';

    my $sub = $filter_hash->{test};

    my($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF);
    $CHROM  = [split("\t",$test1)]->[0];
    $POS    = [split("\t",$test1)]->[1];
    $FAIL   = 1;
    $PASS   = 0;
    #$VCF    = $$opts{vcf};

    is($filter_hash->{tag}, 'INFO/LEN',"_test_FF005 check the correct info tag has been set for the rule");

    $RECORD = [split("\t",$test1)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_FF005 depth < 200");
    $RECORD = [split("\t",$test2)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_FF005 depth < 200 with different pos/neg ratios");
    $RECORD = [split("\t",$test3)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_FF005 depth > 200");
    $RECORD = [split("\t",$test4)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_FF005 depth > 200 with different pos/neg ratios");
    $RECORD = [split("\t",$test5)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_FF005 depth > 200");
    $RECORD = [split("\t",$test6)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_FF005 depth > 200 with different pos/neg ratios");
    $RECORD = [split("\t",$test7)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_FF005 depth > 200 pMtPos == 1 pMtNeg == 1");
    $RECORD = [split("\t",$test8)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_FF005 depth > 200 pMtPos == 4 pMtNeg == 3");
    $RECORD = [split("\t",$test9)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_FF005 depth > 200 pMtPos == 4 pMtNeg == 4");
    $RECORD = [split("\t",$test10)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_FF005 depth > 200 pMtPos == 4 pMtNeg == 5");
    $RECORD = [split("\t",$test11)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_FF005 depthNeg > 200 pMtPos == 1");
    $RECORD = [split("\t",$test12)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_FF005  depthPos > 200 pMtPos == 7");
    $RECORD = [split("\t",$test13)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_FF005  depthPos > 200 pMtPos == 8");
    $RECORD = [split("\t",$test14)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_FF005 depthPos > 200 pMtNeg == 1");
    $RECORD = [split("\t",$test15)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_FF005  depthNeg > 200 pMtNeg == 7");
    $RECORD = [split("\t",$test16)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_FF005  depthNeg > 200 pMtNeg == 8");
	};
}

sub _test_FF006{
	my ($filter_hash) = @_;
  subtest "Test rule FF006" => sub {
    my $test1 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=4;SM=138;S1=10;S2=203.791;REP=9	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:0:0:0:0	0:0:0:0:0:0:0:0:0:0';
    my $test2 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=4;SM=138;S1=10;S2=203.791;REP=10	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:0:0:0:0	0:0:0:0:0:0:0:0:0:0';
    my $test3 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=4;SM=138;S1=10;S2=203.791;REP=0	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:0:0:0:0	0:0:0:0:0:0:0:0:0:0';
    my $test4 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=5;SM=138;S1=10;S2=203.791;REP=9	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:0:0:0:0	0:0:0:0:0:0:0:0:0:0';
    my $test5 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=5;SM=138;S1=10;S2=203.791;REP=10	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:0:0:0:0	0:0:0:0:0:0:0:0:0:0';
    my $test6 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=5;SM=138;S1=10;S2=203.791;REP=0	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:0:0:0:0	0:0:0:0:0:0:0:0:0:0';
    my $test7 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=4;SM=138;S1=10;S2=203.791;REP=11	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:0:0:0:0	0:0:0:0:0:0:0:0:0:0';

    my $sub = $filter_hash->{test};

    my($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF);
    $CHROM  = [split("\t",$test1)]->[0];
    $POS    = [split("\t",$test1)]->[1];
    $FAIL   = 1;
    $PASS   = 0;
    #$VCF    = $$opts{vcf};

    is($filter_hash->{tag}, 'INFO/LEN',"_test_FF006 check the correct info tag has been set for the rule");

    $RECORD = [split("\t",$test1)];
    $MATCH = 4;
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_FF006 length == 4 rep == 9");
    $RECORD = [split("\t",$test2)];
    $MATCH = 4;
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_FF006 length == 4 rep == 10");
    $RECORD = [split("\t",$test3)];
    $MATCH = 4;
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_FF006 length == 4 rep == 0");
    $RECORD = [split("\t",$test7)];
    $MATCH = 4;
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_FF006 length == 4 rep == 11");
    $RECORD = [split("\t",$test4)];
    $MATCH = 5;
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_FF006 length == 5 rep == 9");
    $RECORD = [split("\t",$test5)];
    $MATCH = 5;
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_FF006 length == 5 rep == 10");
    $RECORD = [split("\t",$test6)];
    $MATCH = 5;
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_FF006 length == 5 rep == 0");
  };
}

sub _test_FF007{
	my ($filter_hash) = @_;
  subtest "Test rule FF007" => sub {
    my $test1 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	FD	0	5';
    my $test2 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	FD	1	5';
    my $test3 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	FD	0	6';
    my $test4 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	FD	1	6';
    my $test5 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	FD	7	100';
    my $test6 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	FD	8	100';
    my $test7 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	FD	9	100';

    my $sub = $filter_hash->{test};

    my($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF);
    $CHROM  = [split("\t",$test1)]->[0];
    $POS    = [split("\t",$test1)]->[1];
    $FAIL   = 1;
    $PASS   = 0;
    #$VCF    = $$opts{vcf};

    is($filter_hash->{tag}, 'INFO/LEN',"_test_FF007 check the correct info tag has been set for the rule");

    $RECORD = [split("\t",$test1)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_FF007 rdMt < 6");
    $RECORD = [split("\t",$test2)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_FF007 rdMt < 6 and rdWt > 8pc");

    $RECORD = [split("\t",$test3)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_FF007 rdMt == 6 and rdWt < 8pc");
    $RECORD = [split("\t",$test4)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_FF007 rdMt == 6 and rdWt > 8pc");

    $RECORD = [split("\t",$test5)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_FF007 rdMt == 100 and rdWt == 7");
    $RECORD = [split("\t",$test6)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_FF007 rdMt == 100 and rdWt == 8");
    $RECORD = [split("\t",$test7)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_FF007 rdMt == 100 and rdWt > 8");
	};
}

sub _test_FF008{
	my ($filter_hash) = @_;
  subtest "Test rule FF008" => sub {
    my $test1 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	FC	6	100';
    my $test2 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	FC	5	100';
    my $test3 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	FC	4	100';
    my $test4 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	FC	0	100';

    my $sub = $filter_hash->{test};

    my($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF);
    $CHROM  = [split("\t",$test1)]->[0];
    $POS    = [split("\t",$test1)]->[1];
    $FAIL   = 1;
    $PASS   = 0;
    #$VCF    = $$opts{vcf};

    is($filter_hash->{tag}, 'INFO/REP',"_test_FF008 check the correct info tag has been set for the rule");

    $RECORD = [split("\t",$test1)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_FF008 pWt > pMt5%");
    $RECORD = [split("\t",$test2)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_FF008 pWt == pMt5%");
    $RECORD = [split("\t",$test3)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_FF008 pWt < pMt5% ");
    $RECORD = [split("\t",$test4)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_FF008 pWt == 0 ");
	};
}

sub _test_FF009{
	#fail; #is coding filter
}

sub _test_FF010{
	#fail; #normal panel filter
}

sub _test_FF011{
	#fail; #microsatelite filter
}

sub _test_FF012{
	my ($filter_hash) = @_;
  subtest "Test rule FF0012" => sub {
    my $test1 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	FC:FD:PD:ND	0:0:9:0	0:0:9:0';
    my $test2 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	FC:FD:PD:ND	0:0:0:9	0:0:0:9';
    my $test3 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	FC:FD:PD:ND	0:0:10:0	0:0:10:0';
    my $test4 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	FC:FD:PD:ND	0:0:0:10	0:0:0:10';
    my $test5 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	FC:FD:PD:ND	0:0:9:0	0:0:10:0';
    my $test6 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	FC:FD:PD:ND	0:0:10:0	0:0:9:0';
    my $test7 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	FC:FD:PD:ND	2:10:10:0	2:10:10:0';
    my $test8 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	FC:FD:PD:ND	1:10:10:0	1:10:10:0';
    my $test9 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	FC:FD:PD:ND	3:10:10:0	3:10:10:0';
    my $test10 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	FC:FD:PD:ND	2:10:10:0	1:10:10:0';
    my $test11 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	FC:FD:PD:ND	1:10:10:0	2:10:10:0';

    my $sub = $filter_hash->{test};

    my($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF);
    $MATCH  = 1;
    $CHROM  = [split("\t",$test1)]->[0];
    $POS    = [split("\t",$test1)]->[1];
    $FAIL   = 1;
    $PASS   = 0;
    #$VCF    = $$opts{vcf};

    is($filter_hash->{tag}, 'INFO/LEN',"_test_F012 check the correct info tag has been set for the rule");

    $RECORD = [split("\t",$test1)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F012 dWt == 9 dMt == 9");
    $RECORD = [split("\t",$test2)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F012 dWt == 9 dMt == 9 with different pos/neg ratios");
    $RECORD = [split("\t",$test3)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_F012 dWt == 10 dMt == 10");
    $RECORD = [split("\t",$test4)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_F012 dWt == 10 dMt == 10 with different pos/neg ratios");
    $RECORD = [split("\t",$test5)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F012 dWt == 9 dMt == 10");
    $RECORD = [split("\t",$test6)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F012 dWt == 10 dMt == 9");
    $RECORD = [split("\t",$test7)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_F012 WtFC == 2 WtFD == 10 MtFC == 2 MtFD == 10");
    $RECORD = [split("\t",$test8)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F012 WtFC == 1 WtFD == 10 MtFC == 1 MtFD == 10");
    $RECORD = [split("\t",$test9)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_F012 WtFC == 3 WtFD == 10 MtFC == 3 MtFD == 10");
    $RECORD = [split("\t",$test10)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F012 WtFC == 2 WtFD == 10 MtFC == 1 MtFD == 10");
    $RECORD = [split("\t",$test11)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F012 WtFC == 1 WtFD == 10 MtFC == 2 MtFD == 10");
  }
}

sub _test_FF015{
	my ($filter_hash) = @_;
  subtest "Test rule FF015" => sub {
    my $test1 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	FC	0	0';
    my $test2 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	FC	1	0';

    my $sub = $filter_hash->{test};

    my($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF);
    $CHROM  = [split("\t",$test1)]->[0];
    $POS    = [split("\t",$test1)]->[1];
    $FAIL   = 1;
    $PASS   = 0;
    #$VCF    = $$opts{vcf};

    is($filter_hash->{tag}, 'INFO/LEN',"_test_FF015 check the correct info tag has been set for the rule");
    $RECORD = [split("\t",$test1)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_FF0015 no wild type at all");
    $RECORD = [split("\t",$test2)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_FF0015 pWtFC == 1");
  };
}

sub _test_FF016{
	my ($filter_hash) = @_;
  subtest "Test rule FF016" => sub {
    my $test1 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=1	PP:NP:PB:NB	0:0:0:0	2:3:1:0';
    my $test2 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=1	PP:NP:PB:NB	0:0:0:0	5:0:1:0';
    my $test3 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=1	PP:NP:PB:NB	0:0:0:0	0:5:1:0';
    my $test4 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=1	PP:NP:PB:NB	0:0:0:0	6:0:1:0';
    my $test5 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=1	PP:NP:PB:NB	0:0:0:0	4:0:0:0';
    my $test6 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=1	PP:NP:PB:NB	0:0:0:0	0:4:0:0';
    my $test7 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=1	PP:NP:PB:NB	0:0:0:0	5:1:1:0';
     my $test8 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=1	PP:NP:PB:NB	0:0:0:0	0:5:0:0';
     my $test9 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=1	PP:NP:PB:NB	0:0:0:0	4:1:0:0';
    my $test10 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=1	PP:NP:PB:NB	0:0:0:0	1:4:0:0';
    my $test11 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=2	PP:NP:PB:NB	0:0:0:0	0:5:0:0';
    my $test12 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=2	PP:NP:PB:NB	0:0:0:0	4:1:0:0';
    my $test13 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=2	PP:NP:PB:NB	0:0:0:0	1:4:0:0';

    my $test14 = '22	16404839	.	G	GA	.	.	PC=I;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=0	PP:NP:PB:NB	0:0:0:0	1:4:0:0';
    my $test15 = '22	16404839	.	G	GA	.	.	PC=I;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=1	PP:NP:PB:NB	0:0:0:0	1:4:0:0';


    my $sub = $filter_hash->{test};

    my($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF);
    $CHROM  = [split("\t",$test1)]->[0];
    $POS    = [split("\t",$test1)]->[1];
    $FAIL   = 1;
    $PASS   = 0;
    #$VCF    = $$opts{vcf};

    is($filter_hash->{tag}, 'INFO/REP',"_test_FF016 check the correct info tag has been set for the rule");

    $RECORD = [split("\t",$test1)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_FF016 pMt == 5");
    $RECORD = [split("\t",$test2)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_FF016 pMt == 5 with different pos/neg ratios");
    $RECORD = [split("\t",$test3)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_FF016 pMt == 5 with different pos/neg ratios");
    $RECORD = [split("\t",$test4)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_FF016 pMt == 6");
    $RECORD = [split("\t",$test5)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_FF016 pMt == 4");
    $RECORD = [split("\t",$test6)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_FF016 pMt == 4 with different pos/neg ratios");
    $RECORD = [split("\t",$test7)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_FF016 pMt == 5 bMt > 0");

    $MATCH = 1;
    $RECORD = [split("\t",$test8)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_FF016 pMt == 5 bMt > 0 rep == 1");
    $MATCH = 1;
    $RECORD = [split("\t",$test9)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_FF016 pMt == 5 bMt > 0 rep == 1 pMtPos > 1 pMtPos == 0");
    $MATCH = 1;
    $RECORD = [split("\t",$test10)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_FF016 pMt == 5 bMt > 0 rep == 1 pMtPos == 0 pMtNeg > 1");

    $MATCH = 2;
    $RECORD = [split("\t",$test11)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_FF016 pMt == 5 bMt > 0 rep == 2: pindel calls in only one direction");
    $MATCH = 2;
    $RECORD = [split("\t",$test12)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_FF016 pMt == 5 bMt > 0 rep == 2: pindel calls in both directions");
    $MATCH = 2;
    $RECORD = [split("\t",$test13)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_FF016 pMt == 5 bMt > 0 rep == 2: pindel calls in both direction");

    $MATCH = 0;
    $RECORD = [split("\t",$test14)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_FF016 pMt == 5 bMt > 0 rep == 0: insertion pindel calls in both directions");
    $MATCH = 1;
    $RECORD = [split("\t",$test15)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_FF016 pMt == 5 bMt > 0 rep == 1: insertion pindel calls in both directions");
  };
}

sub _test_FF017{
	#fail; #simple repeat filter
}

sub _test_FF018{
	my ($filter_hash) = @_;
  subtest "Test rule FF018" => sub {
    my $test1 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PR:NR	5:5	5:5';
    my $test2 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PR:NR	6:4	5:5';
    my $test3 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PR:NR	5:5	6:4';
    my $test4 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PR:NR	5:4	5:5';
    my $test5 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PR:NR	5:5	5:4';
    my $test6 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PR:NR	5:6	5:5';
    my $test7 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PR:NR	5:5	5:6';

    my $sub = $filter_hash->{test};

    my($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF);
    $CHROM  = [split("\t",$test1)]->[0];
    $POS    = [split("\t",$test1)]->[1];
    $FAIL   = 1;
    $PASS   = 0;
    #$VCF    = $$opts{vcf};

    is($filter_hash->{tag}, 'INFO/LEN',"_test_FF018 check the correct info tag has been set for the rule");

    $RECORD = [split("\t",$test1)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_FF014 pWt == 10 pMt == 10");
    $RECORD = [split("\t",$test2)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_FF014 pWt == 10 pMt == 10 with different pos/neg ratios");
    $RECORD = [split("\t",$test3)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_FF014 pWt == 10 pMt == 10 with different pos/neg ratios");
    $RECORD = [split("\t",$test4)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_FF014 pWt == 10 pMt == 9");
    $RECORD = [split("\t",$test5)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_FF014 pWt == 9 pMt == 10");
    $RECORD = [split("\t",$test6)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_FF014 pWt == 11 pMt == 10");
    $RECORD = [split("\t",$test7)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_FF014 pWt == 10 pMt == 11");
  };
}

sub _test_FF019{
	my ($filter_hash) = @_;
  subtest "Test rule FF019" => sub {
    my $test1 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	FC:FD	0:0	5:20';
    my $test2 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	FC:FD	0:0	3:60';
    my $test3 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	FC:FD	0:0	3:59';
    my $test4 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	FC:FD	0:0	1:20';
    my $test5 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	FC:FD	0:0	3:100';
    my $test6 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	FC:FD	3:59	3:60';

    my $sub = $filter_hash->{test};

    my($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF);
    $CHROM  = [split("\t",$test1)]->[0];
    $POS    = [split("\t",$test1)]->[1];
    $FAIL   = 1;
    $PASS   = 0;
    #$VCF    = $$opts{vcf};

    is($filter_hash->{tag}, 'INFO/LEN',"_test_FF019 check the correct info tag has been set for the rule");

    $RECORD = [split("\t",$test1)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_FF019 pWtFC == 0 pMtFC == 5 pMtFD == 20 pMtFC/pMtFD > 0.05");
    $RECORD = [split("\t",$test2)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_FF019 pWtFC == 0 pMtFC == 3 pMtFD == 60 pMtFC/pMtFD == 0.05");
    $RECORD = [split("\t",$test3)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_FF019 pWtFC == 0 pMtFC == 3 pMtFD == 59 pMtFC/pMtFD > 0.05");
    $RECORD = [split("\t",$test4)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_FF019 pWtFC == 0 pMtFC == 1 pMtFD == 20 pMtFC/pMtFD == 0.05");
    $RECORD = [split("\t",$test5)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_FF019 pWtFC == 0 pMtFC == 3 pMtFD == 61 pMtFC/pMtFC < 0.05");
    $RECORD = [split("\t",$test6)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_FF019 pWtFC / PWtFD < 0.05 pMtFC / pMtFD == 0.05");
  };
}

sub _test_FF020{
	my ($filter_hash) = @_;
  subtest "Test rule FF020" => sub {
    my $test1 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	FC:FD	5:200	19:100';
    my $test2 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	FC:FD	3:200	19:100';
    my $test3 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	FC:FD	5:200	21:100';
    my $test4 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	FC:FD	1:10	20:99';
    my $test5 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	FC:FD	0:10	20:101';
    my $test6 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	FC:FD	1:11	20:111';
    my $test7 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	FC:FD	1:9	20:91';
    my $test8 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	FC:FD	2:40	20:100';
    my $test9 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	FC:FD	2:40	19:100';
    my $test10 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	FC:FD	2:41	20:100';
    my $test11 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	FC:FD	2:41	20:99';
    my $test12 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	FC:FD	2:40	20:99';
    my $test13 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	FC:FD	2:39	20:99';
    my $test14 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	FC:FD	1:11	9:111';
    my $test15 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	FC:FD	2:100	20:101';
    my $test16 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	FC:FD	3:10	10:100';
    my $test17 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	FC:FD	0:0	0:0';
    my $test18 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	FC:FD	4:200	1:100';
    my $test19 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	FC:FD	3:200	20:100';
    my $test20 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	FC:FD	5:201	19:100';

    my $sub = $filter_hash->{test};

    my($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF);
    $CHROM  = [split("\t",$test1)]->[0];
    $POS    = [split("\t",$test1)]->[1];
    $FAIL   = 1;
    $PASS   = 0;
    #$VCF    = $$opts{vcf};

    is($filter_hash->{tag}, 'INFO/LEN',"_test_FF019 check the correct info tag has been set for the rule");

    $RECORD = [split("\t",$test1)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_FF020 WtFD == 200 WtFC / WtFD > 0.02 tumFC / tumFD < 0.2");
    $RECORD = [split("\t",$test2)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_FF020 WtFD == 200 WtFC / WtFD < 0.02 tumFC / tumFD < 0.2");
    $RECORD = [split("\t",$test3)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_FF020 WtFD == 200 WtFC / WtFD > 0.02 tumFC / tumFD > 0.2");
    $RECORD = [split("\t",$test4)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_FF020 WtFD < 200 WtFC == 1 WtFD == 10 WtFC < MtFC * 0.1);");
    $RECORD = [split("\t",$test5)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_FF020 WtFD < 200 WtFC == 0 WtFD == 10 WtFC < MtFC * 0.1);");
    $RECORD = [split("\t",$test6)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_FF020 WtFD < 200 WtFC == 1 WtFD == 11 WtFC < MtFC * 0.1);");
    $RECORD = [split("\t",$test7)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_FF020 WtFD < 200 WtFC == 1 WtFD == 9 WtFC / WtFD > 0.05);");
    $RECORD = [split("\t",$test8)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_FF020 WtFD < 200 WtFC == 2 WtFC / WtFD == 0.05 MtFC / MtFD == 0.2);");
    $RECORD = [split("\t",$test9)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_FF020 WtFD < 200 WtFC == 2 WtFC / WtFD == 0.05 MtFC / MtFD < 0.2);");
    $RECORD = [split("\t",$test10)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_FF020 WtFD < 200 WtFC == 2 WtFC / WtFD < 0.05 MtFC / MtFD == 0.2);");
    $RECORD = [split("\t",$test11)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_FF020 WtFD < 200 WtFC == 2 WtFC / WtFD < 0.05 MtFC / MtFD > 0.2);");
    $RECORD = [split("\t",$test12)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_FF020 WtFD < 200 WtFC == 2 WtFC / WtFD == 0.05 MtFC / MtFD > 0.2);");
    $RECORD = [split("\t",$test13)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_FF020 WtFD < 200 WtFC == 2 WtFC / WtFD > 0.05 MtFC / MtFD > 0.2);");
    $RECORD = [split("\t",$test14)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_FF020 WtFD < 200 WtFC == 1 WtFD == 11 WtFC > MtFC * 0.1 MtFC / MtFD < 0.2);");
    $RECORD = [split("\t",$test15)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_FF020 WtFD < 200 WtFC == 1 WtFD == 100 WtFC == MtFC * 0.1 MtFC / MtFD < 0.2);");
    $RECORD = [split("\t",$test16)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_FF020 WtFD < 200 WtFC == 3");
    $RECORD = [split("\t",$test17)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_FF020 divid by 0");
    $RECORD = [split("\t",$test18)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_FF020 WtFD == 200 WtFC / WtFD == 0.02 tumFC / tumFD < 0.2");
    $RECORD = [split("\t",$test19)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_FF020 WtFD == 200 WtFC / WtFD < 0.02 tumFC / tumFD == 0.2");
    $RECORD = [split("\t",$test20)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_FF020 WtFD > 200 WtFC / WtFD > 0.02 tumFC / tumFD < 0.2");
  };
}
