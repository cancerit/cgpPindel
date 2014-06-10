####################################################
# Copyright (c) 2013 Genome Research Ltd.
# Author: Cancer Genome Project, cgpit@sanger.ac.uk
# See LICENCE.TXT for details
####################################################

use strict;
use warnings;
use Cwd 'abs_path';
use Test::More;
use Data::Dumper;
use Const::Fast qw(const);

const my @AVAILABLE_RULES => qw(F001 F002 F003 F004 F005 F006 F007 F008 F009 F010 F012 F015 F016 F017 F018);

my %rule_test_dispatch = ('F001' => \&_test_F001,
                          'F002' => \&_test_F002,
                          'F003' => \&_test_F003,
                          'F004' => \&_test_F004,
                          'F005' => \&_test_F005,
                          'F006' => \&_test_F006,
                          'F007' => \&_test_F007,
                          'F008' => \&_test_F008,
                          'F009' => \&_test_F009,
                          'F010' => \&_test_F010,
                          'F012' => \&_test_F012,
                          'F015' => \&_test_F015,
                          'F016' => \&_test_F016,
                          'F017' => \&_test_F017,
                          'F018' => \&_test_F018,
                        );

use_ok('Sanger::CGP::PindelPostProcessing::VcfSoftFlagger');
use_ok('Sanger::CGP::PindelPostProcessing::FilterRules');

my @rules_found = Sanger::CGP::PindelPostProcessing::FilterRules::available_rules();
is_deeply(\@rules_found, \@AVAILABLE_RULES, 'Expected set of rules are implemented');
for my $flag(@AVAILABLE_RULES) {
  $rule_test_dispatch{$flag}(Sanger::CGP::PindelPostProcessing::FilterRules->rule($flag));
}

done_testing();

sub _test_F001{
	my ($filter_hash) = @_;
  subtest "Test rule F001" => sub {
    my $test1 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	5:5:0:0:0:0:0:0:0:0	5:5:0:0:0:0:0:0:0:0';
    my $test2 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	4:6:0:0:0:0:0:0:0:0	6:4:0:0:0:0:0:0:0:0';
    my $test3 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	5:6:0:0:0:0:0:0:0:0	5:5:0:0:0:0:0:0:0:0';
    my $test4 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	5:4:0:0:0:0:0:0:0:0	5:5:0:0:0:0:0:0:0:0';

    my $sub = $filter_hash->{test};

    my($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF);
    $CHROM  = [split("\t",$test1)]->[0];
    $POS    = [split("\t",$test1)]->[1];
    $FAIL   = 1;
    $PASS   = 0;
    #$VCF    = $$opts{vcf};

    is($filter_hash->{tag}, 'INFO/LEN',"_test_F001 check the correct info tag has been set for the rule");

    $RECORD = [split("\t",$test1)];

    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_F001 pWt == pMt");
    $RECORD = [split("\t",$test2)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_F001 pWt == pMt with different pos/neg ratios");
    $RECORD = [split("\t",$test3)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_F001 pWt > pMt");
    $RECORD = [split("\t",$test4)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F001 pWt < pMt");
	};
}


sub _test_F002{
	my ($filter_hash) = @_;
  subtest "Test rule F002" => sub {
    my $test1 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=4;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	1:0:0:0:0:0:0:0:0:0	0:0:0:0:0:0:0:0:0:0';
    my $test2 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=4;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:1:0:0:0:0:0:0:0:0	0:0:0:0:0:0:0:0:0:0';
    my $test3 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=4;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:0:0:0:0	0:0:0:0:0:0:0:0:0:0';

    my $test4 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=5;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	1:0:0:0:0:0:0:0:0:0	0:0:0:0:0:0:0:0:0:0';
    my $test5 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=5;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:1:0:0:0:0:0:0:0:0	0:0:0:0:0:0:0:0:0:0';
    my $test6 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=5;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:0:0:0:0	0:0:0:0:0:0:0:0:0:0';


    my $sub = $filter_hash->{test};

    my($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF);
    $CHROM  = [split("\t",$test1)]->[0];
    $POS    = [split("\t",$test1)]->[1];
    $FAIL   = 1;
    $PASS   = 0;
    #$VCF    = $$opts{vcf};

    is($filter_hash->{tag}, 'INFO/LEN',"_test_F002 check the correct info tag has been set for the rule");

    $RECORD = [split("\t",$test1)];
    $MATCH = 4;
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F002 length == 4 pWt == 1");
    $RECORD = [split("\t",$test2)];
    $MATCH = 4;
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F002 length == 4 pWt == 1 different pos/neg ratios");
    $RECORD = [split("\t",$test3)];
    $MATCH = 4;
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F002 length == 4 pWt == 0");

    $RECORD = [split("\t",$test4)];
    $MATCH = 5;
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_F002 length == 5 pWt == 1");
    $RECORD = [split("\t",$test5)];
    $MATCH = 5;
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_F002 length == 5 pWt == 1 different pos/neg ratios");
    $RECORD = [split("\t",$test6)];
    $MATCH = 5;
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F002 length == 5 pWt == 0");
  };
}

sub _test_F003{
	my ($filter_hash) = @_;
  subtest "Test rule F003" => sub {
    my $test1 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:0:0:0:0	3:0:0:0:0:0:0:0:0:0';
    my $test2 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:0:0:0:0	0:3:0:0:0:0:0:0:0:0';
    my $test3 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:0:0:0:0	4:0:0:0:0:0:0:0:0:0';
    my $test4 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:0:0:0:0	0:4:0:0:0:0:0:0:0:0';
    my $test5 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:0:0:0:0	3:4:0:0:0:0:0:0:0:0';
    my $test6 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:0:0:0:0	2:2:0:0:0:0:0:0:0:0';
    my $test7 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:0:0:0:0	2:1:0:0:0:0:0:0:0:0';
    my $test8 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:0:0:0:0	1:2:0:0:0:0:0:0:0:0';
    my $test9 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:0:0:0:0	2:0:0:0:0:0:0:0:0:0';
    my $test10 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:0:0:0:0	0:2:0:0:0:0:0:0:0:0';
    my $test11 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:0:0:0:0	0:0:0:0:0:0:0:0:0:0';



    my $sub = $filter_hash->{test};

    my($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF);
    $CHROM  = [split("\t",$test1)]->[0];
    $POS    = [split("\t",$test1)]->[1];
    $FAIL   = 1;
    $PASS   = 0;
    #$VCF    = $$opts{vcf};

    is($filter_hash->{tag}, 'INFO/LEN',"_test_F003 check the correct info tag has been set for the rule");

    $RECORD = [split("\t",$test1)];

    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F003 pMtPos == 3 pMtNeg == 0");
    $RECORD = [split("\t",$test2)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F003 pMtPos == 0 pMtNeg == 3");
    $RECORD = [split("\t",$test3)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F003 pMtPos == 4 pMtNeg == 0");
    $RECORD = [split("\t",$test4)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F003 pMtPos == 0 pMtNeg == 4");

    $RECORD = [split("\t",$test5)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F003 pMtPos == 3 pMtNeg == 4");
    $RECORD = [split("\t",$test6)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F003 pMtPos == 2 pMtNeg == 2");
    $RECORD = [split("\t",$test7)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_F003 pMtPos == 2 pMtNeg == 1");
    $RECORD = [split("\t",$test8)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_F003 pMtPos == 1 pMtNeg == 2");
    $RECORD = [split("\t",$test9)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_F003 pMtPos == 2 pMtNeg == 0");
    $RECORD = [split("\t",$test10)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_F003 pMtPos == 0 pMtNeg == 2");
    $RECORD = [split("\t",$test11)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_F003 pMtPos == 0 pMtNeg == 0");
	};
}


sub _test_F004{
	my ($filter_hash) = @_;
  subtest "Test rule F004" => sub {
    my $test1 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:0:0:0:0	0:0:0:0:0:0:9:0:0:0';
    my $test2 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:0:0:0:0	0:0:0:0:0:0:10:0:0:0';
    my $test3 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:0:0:0:0	0:0:0:0:0:0:199:0:0:0';
    my $test4 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:0:0:0:0	0:0:0:0:0:0:200:0:0:0';

    my $test5 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:0:0:0:0	0:0:0:0:0:0:100:0:2:3';
    my $test6 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:0:0:0:0	0:0:0:0:0:0:100:0:2:2';


    my $test7 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:0:0:0:0	0:0:0:0:0:0:0:100:4:0';
    my $test8 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:0:0:0:0	0:0:0:0:0:0:100:0:8:0';
    my $test9 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:0:0:0:0	0:0:0:0:0:0:100:0:7:0';


    my $test10 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:0:0:0:0	0:0:0:0:0:0:100:0:0:4';
    my $test11 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:0:0:0:0	0:0:0:0:0:0:0:100:0:8';
    my $test12 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:0:0:0:0	0:0:0:0:0:0:0:100:0:7';


    my $sub = $filter_hash->{test};

    my($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF);
    $CHROM  = [split("\t",$test1)]->[0];
    $POS    = [split("\t",$test1)]->[1];
    $FAIL   = 1;
    $PASS   = 0;
    #$VCF    = $$opts{vcf};

    is($filter_hash->{tag}, 'INFO/LEN',"_test_F004 check the correct info tag has been set for the rule");

    $RECORD = [split("\t",$test1)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F004 ");
    $RECORD = [split("\t",$test2)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_F004 ");
    $RECORD = [split("\t",$test3)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_F004 ");
    $RECORD = [split("\t",$test4)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F004 ");

    $RECORD = [split("\t",$test5)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F004 ");
    $RECORD = [split("\t",$test6)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_F004 ");

    $RECORD = [split("\t",$test7)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F004 ");

    $RECORD = [split("\t",$test8)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F004 ");
    $RECORD = [split("\t",$test9)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_F004 ");

    $RECORD = [split("\t",$test10)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F004 ");

    $RECORD = [split("\t",$test11)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F004 ");
    $RECORD = [split("\t",$test12)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_F004 ");
  };
}

sub _test_F005{
	my ($filter_hash) = @_;
  subtest "Test rule F005" => sub {
    my $test1 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:0:0:0:0	0:0:0:0:0:0:199:0:1:0';
    my $test2 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:0:0:0:0	0:0:0:0:0:0:0:199:1:0';

    my $test3 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:0:0:0:0	0:0:0:0:0:0:200:0:4:4';
    my $test4 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:0:0:0:0	0:0:0:0:0:0:0:200:4:4';
    my $test5 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:0:0:0:0	0:0:0:0:0:0:199:1:4:4';
    my $test6 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:0:0:0:0	0:0:0:0:0:0:1:199:4:4';

    my $test7 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:0:0:0:0	0:0:0:0:0:0:200:0:1:1';
    my $test8 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:0:0:0:0	0:0:0:0:0:0:200:0:3:4';
    my $test9 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:0:0:0:0	0:0:0:0:0:0:200:0:4:4';
    my $test10 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:0:0:0:0	0:0:0:0:0:0:200:0:4:5';

    my $test11 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:0:0:0:0	0:0:0:0:0:0:0:200:1:0';
    my $test12 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:0:0:0:0	0:0:0:0:0:0:200:0:7:0';
    my $test13 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:0:0:0:0	0:0:0:0:0:0:200:0:8:0';

    my $test14 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:0:0:0:0	0:0:0:0:0:0:200:0:0:1';
    my $test15 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:0:0:0:0	0:0:0:0:0:0:0:200:0:7';
    my $test16 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:0:0:0:0	0:0:0:0:0:0:0:200:0:8';

    my $sub = $filter_hash->{test};

    my($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF);
    $CHROM  = [split("\t",$test1)]->[0];
    $POS    = [split("\t",$test1)]->[1];
    $FAIL   = 1;
    $PASS   = 0;
    #$VCF    = $$opts{vcf};

    is($filter_hash->{tag}, 'INFO/LEN',"_test_F005 check the correct info tag has been set for the rule");

    $RECORD = [split("\t",$test1)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F005 depth < 200");
    $RECORD = [split("\t",$test2)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F005 depth < 200 with different pos/neg ratios");
    $RECORD = [split("\t",$test3)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F005 depth > 200");
    $RECORD = [split("\t",$test4)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F005 depth > 200 with different pos/neg ratios");
    $RECORD = [split("\t",$test5)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F005 depth > 200");
    $RECORD = [split("\t",$test6)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F005 depth > 200 with different pos/neg ratios");
    $RECORD = [split("\t",$test7)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_F005 depth > 200 pMtPos == 1 pMtNeg == 1");
    $RECORD = [split("\t",$test8)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_F005 depth > 200 pMtPos == 4 pMtNeg == 3");
    $RECORD = [split("\t",$test9)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F005 depth > 200 pMtPos == 4 pMtNeg == 4");
    $RECORD = [split("\t",$test10)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F005 depth > 200 pMtPos == 4 pMtNeg == 5");
    $RECORD = [split("\t",$test11)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F005 depthNeg > 200 pMtPos == 1");
    $RECORD = [split("\t",$test12)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_F005  depthPos > 200 pMtPos == 7");
    $RECORD = [split("\t",$test13)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F005  depthPos > 200 pMtPos == 8");
    $RECORD = [split("\t",$test14)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F005 depthPos > 200 pMtNeg == 1");
    $RECORD = [split("\t",$test15)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_F005  depthNeg > 200 pMtNeg == 7");
    $RECORD = [split("\t",$test16)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F005  depthNeg > 200 pMtNeg == 8");
	};
}

sub _test_F006{
	my ($filter_hash) = @_;
  subtest "Test rule F006" => sub {
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

    is($filter_hash->{tag}, 'INFO/LEN',"_test_F006 check the correct info tag has been set for the rule");

    $RECORD = [split("\t",$test1)];
    $MATCH = 4;
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F006 length == 4 rep == 9");
    $RECORD = [split("\t",$test2)];
    $MATCH = 4;
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_F006 length == 4 rep == 10");
    $RECORD = [split("\t",$test3)];
    $MATCH = 4;
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F006 length == 4 rep == 0");
    $RECORD = [split("\t",$test7)];
    $MATCH = 4;
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_F006 length == 4 rep == 11");
    $RECORD = [split("\t",$test4)];
    $MATCH = 5;
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F006 length == 5 rep == 9");
    $RECORD = [split("\t",$test5)];
    $MATCH = 5;
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F006 length == 5 rep == 10");
    $RECORD = [split("\t",$test6)];
    $MATCH = 5;
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F006 length == 5 rep == 0");
  };
}

sub _test_F007{
	my ($filter_hash) = @_;
  subtest "Test rule F007" => sub {
    my $test1 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:0:0:0:0	0:0:0:0:0:0:5:0:0:0';
    my $test2 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:1:0:0:0	0:0:0:0:0:0:5:0:0:0';
    my $test3 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:0:1:0:0	0:0:0:0:0:0:5:0:0:0';
    my $test4 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:0:0:0:0	0:0:0:0:0:0:6:0:0:0';
    my $test5 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:1:0:0:0	0:0:0:0:0:0:6:0:0:0';
    my $test6 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:0:1:0:0	0:0:0:0:0:0:6:0:0:0';
    my $test7 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:0:0:0:0	0:0:0:0:0:0:3:3:0:0';
    my $test8 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:1:0:0:0	0:0:0:0:0:0:3:3:0:0';
    my $test9 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:0:1:0:0	0:0:0:0:0:0:3:3:0:0';
    my $test10 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:7:0:0:0	0:0:0:0:0:0:50:50:0:0';
    my $test11 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:3:4:0:0	0:0:0:0:0:0:50:50:0:0';
    my $test12 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:8:0:0:0	0:0:0:0:0:0:50:50:0:0';
    my $test13 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:4:4:0:0	0:0:0:0:0:0:50:50:0:0';
    my $test14 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:5:4:0:0	0:0:0:0:0:0:50:50:0:0';

    my $sub = $filter_hash->{test};

    my($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF);
    $CHROM  = [split("\t",$test1)]->[0];
    $POS    = [split("\t",$test1)]->[1];
    $FAIL   = 1;
    $PASS   = 0;
    #$VCF    = $$opts{vcf};

    is($filter_hash->{tag}, 'INFO/LEN',"_test_F007 check the correct info tag has been set for the rule");

    $RECORD = [split("\t",$test1)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_F007 rdMt < 6");
    $RECORD = [split("\t",$test2)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_F007 rdMt < 6 and rdWt > 8pc");
    $RECORD = [split("\t",$test3)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_F007 rdMt < 6 and rdWt > 8pc with different pos/neg ratios");
    $RECORD = [split("\t",$test4)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_F007 rdMt == 6 and rdWt < 8pc");
    $RECORD = [split("\t",$test5)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F007 rdMt == 6 and rdWt > 8pc");
    $RECORD = [split("\t",$test6)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F007 rdMt == 6 and rdWt > 8pc with different pos/neg ratios");
    $RECORD = [split("\t",$test7)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_F007 rdMt == 6 and rdWt < 8pc");
    $RECORD = [split("\t",$test8)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F007 rdMt == 6 and rdWt > 8pc");
    $RECORD = [split("\t",$test9)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F007 rdMt == 6 and rdWt > 8pc with different pos/neg ratios");
    $RECORD = [split("\t",$test10)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_F007 rdMt == 100 and rdWt == 7");
    $RECORD = [split("\t",$test11)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_F007 rdMt == 100 and rdWt == 7 with different pos/neg ratios");
    $RECORD = [split("\t",$test12)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F007 rdMt == 100 and rdWt == 8");
    $RECORD = [split("\t",$test13)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F007 rdMt == 100 and rdWt == 8 with different pos/neg ratios");
    $RECORD = [split("\t",$test14)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F007 rdMt == 100 and rdWt > 8");
	};
}

sub _test_F008{
	my ($filter_hash) = @_;
  subtest "Test rule F008" => sub {
    my $test1 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	6:0:0:0:0:0:0:0:0:0	50:50:0:0:0:0:0:0:0:0';
    my $test2 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:6:0:0:0:0:0:0:0:0	50:50:0:0:0:0:0:0:0:0';
    my $test3 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	3:3:0:0:0:0:0:0:0:0	50:50:0:0:0:0:0:0:0:0';
    my $test4 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	5:1:0:0:0:0:0:0:0:0	50:50:0:0:0:0:0:0:0:0';
    my $test5 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	5:0:0:0:0:0:0:0:0:0	50:50:0:0:0:0:0:0:0:0';
    my $test6 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:5:0:0:0:0:0:0:0:0	50:50:0:0:0:0:0:0:0:0';
    my $test7 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	1:4:0:0:0:0:0:0:0:0	50:50:0:0:0:0:0:0:0:0';
    my $test8 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	1:3:0:0:0:0:0:0:0:0	50:50:0:0:0:0:0:0:0:0';
    my $test9 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:0:0:0:0	50:50:0:0:0:0:0:0:0:0';

    my $sub = $filter_hash->{test};

    my($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF);
    $CHROM  = [split("\t",$test1)]->[0];
    $POS    = [split("\t",$test1)]->[1];
    $FAIL   = 1;
    $PASS   = 0;
    #$VCF    = $$opts{vcf};

    is($filter_hash->{tag}, 'INFO/REP',"_test_F008 check the correct info tag has been set for the rule");

    $RECORD = [split("\t",$test1)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_F008 pWt > pMt5%");
    $RECORD = [split("\t",$test2)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_F008 pWt > pMt5% with different pos/neg ratios");
    $RECORD = [split("\t",$test3)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_F008 pWt > pMt5% with different pos/neg ratios");
    $RECORD = [split("\t",$test4)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_F008 pWt > pMt5% with different pos/neg ratios");
    $RECORD = [split("\t",$test5)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F008 pWt == pMt5%");
    $RECORD = [split("\t",$test6)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F008 pWt == pMt5% with different pos/neg ratios");
    $RECORD = [split("\t",$test7)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F008 pWt == pMt5% with different pos/neg ratios");
    $RECORD = [split("\t",$test8)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F008 pWt < pMt5% ");
    $RECORD = [split("\t",$test8)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F008 pWt == 0 ");
	};
}

sub _test_F009{
	#fail; #is coding filter
}

sub _test_F010{
	#fail; #normal panel filter
}

sub _test_F011{
	#fail; #microsatelite filter
}

sub _test_F012{
	my ($filter_hash) = @_;
  subtest "Test rule F0012" => sub {
    my $test1 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	2:0:0:0:9:0:0:0:0:0	2:0:0:0:9:0:0:0:0:0';
    my $test2 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	2:0:0:0:0:9:0:0:0:0	2:0:0:0:0:9:0:0:0:0';
    my $test3 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	2:0:0:0:0:10:0:0:0:0	2:0:0:0:10:0:0:0:0:0';
    my $test4 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	2:0:0:0:10:0:0:0:0:0	2:0:0:0:0:10:0:0:0:0';
    my $test5 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	3:0:0:0:0:11:0:0:0:0	3:0:0:0:11:0:0:0:0:0';
    my $test6 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	2:0:0:0:0:10:0:0:0:0	2:0:0:0:9:0:0:0:0:0';
    my $test7 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	20:0:0:0:50:50:0:0:0:0	20:0:0:0:50:50:0:0:0:0';
    my $test8 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	10:10:0:0:50:50:0:0:0:0	10:10:0:0:50:50:0:0:0:0';
    my $test9 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	19:0:0:0:50:50:0:0:0:0	19:0:0:0:50:50:0:0:0:0';
    my $test10 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	18:1:0:0:50:50:0:0:0:0	18:1:0:0:50:50:0:0:0:0';
    my $test11 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	20:0:0:0:50:50:0:0:0:0	19:0:0:0:50:50:0:0:0:0';
    my $test12 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	19:0:0:0:50:50:0:0:0:0	20:0:0:0:50:50:0:0:0:0';
    my $test13 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	21:0:0:0:50:50:0:0:0:0	21:0:0:0:50:50:0:0:0:0';
    my $test14 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:20:0:50:50:0:0:0:0	0:0:20:0:50:50:0:0:0:0';
    my $test15 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:10:10:50:50:0:0:0:0	0:0:10:10:50:50:0:0:0:0';
    my $test16 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:19:0:50:50:0:0:0:0	0:0:19:0:50:50:0:0:0:0';
    my $test17 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:18:1:50:50:0:0:0:0	0:0:18:1:50:50:0:0:0:0';
    my $test18 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:20:0:50:50:0:0:0:0	0:0:19:0:50:50:0:0:0:0';
    my $test19 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:19:0:50:50:0:0:0:0	0:0:20:0:50:50:0:0:0:0';
    my $test20 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:21:0:50:50:0:0:0:0	0:0:21:0:50:50:0:0:0:0';

    my $sub = $filter_hash->{test};

    my($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF);
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
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_F012 dWt > 10 dMt == 10");
    $RECORD = [split("\t",$test6)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F012 dWt > 10 dMt == 9");
    $RECORD = [split("\t",$test7)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_F012 dWt == 100 pWt == 20");
    $RECORD = [split("\t",$test8)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_F012 dWt == 100 pWt == 20 with different pos/neg ratios");
    $RECORD = [split("\t",$test9)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F012 dWt == 100 pWt == 19");
    $RECORD = [split("\t",$test10)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F012 dWt == 100 pWt == 19 with different pos/neg ratios");
    $RECORD = [split("\t",$test11)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F012 dWt == 100 pWt == 20 pMt == 19");
    $RECORD = [split("\t",$test12)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F012 dWt == 100 pWt == 19 pMt == 20");
    $RECORD = [split("\t",$test13)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_F012 dWt == 100 pWt == 21 pMt == 21");
    $RECORD = [split("\t",$test14)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_F012 dWt == 100 bWt == 20");
    $RECORD = [split("\t",$test15)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_F012 dWt == 100 bWt == 20 with different pos/neg ratios");
    $RECORD = [split("\t",$test16)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F012 dWt == 100 bWt == 19");
    $RECORD = [split("\t",$test17)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F012 dWt == 100 bWt == 19 with different pos/neg ratios");
    $RECORD = [split("\t",$test18)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F012 dWt == 100 bWt == 20 bMt == 19");
    $RECORD = [split("\t",$test19)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F012 dWt == 100 bWt == 19 bMt == 20");
    $RECORD = [split("\t",$test20)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_F012 dWt == 100 bWt == 21 bMt == 21");
	};
}

sub _test_F015{
	my ($filter_hash) = @_;
  subtest "Test rule F015" => sub {
    my $test1 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:0:0:0:0	0:0:0:0:0:0:0:0:0:0';
    my $test2 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:0:0:1:0	0:0:0:0:0:0:0:0:0:0';
    my $test3 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:0:0:0:1	0:0:0:0:0:0:0:0:0:0';
    my $test4 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:0:0:1:1	0:0:0:0:0:0:0:0:0:0';

    my $sub = $filter_hash->{test};

    my($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF);
    $CHROM  = [split("\t",$test1)]->[0];
    $POS    = [split("\t",$test1)]->[1];
    $FAIL   = 1;
    $PASS   = 0;
    #$VCF    = $$opts{vcf};

    is($filter_hash->{tag}, 'INFO/LEN',"_test_F015 check the correct info tag has been set for the rule");
    $RECORD = [split("\t",$test1)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F0015 no wild type at all");
    $RECORD = [split("\t",$test2)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_F0015 pWtPos == 1");
    $RECORD = [split("\t",$test3)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_F0015 pWtNeg == 1");
    $RECORD = [split("\t",$test4)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_F0015 pWtPos == 1 bWtPos == 1");
  };
}

sub _test_F016{
	my ($filter_hash) = @_;
  subtest "Test rule F016" => sub {
    my $test1 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=1	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:0:0:0:0	2:3:1:0:0:0:0:0:0:0';
    my $test2 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=1	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:0:0:0:0	5:0:1:0:0:0:0:0:0:0';
    my $test3 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=1	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:0:0:0:0	0:5:1:0:0:0:0:0:0:0';
    my $test4 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=1	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:0:0:0:0	6:0:1:0:0:0:0:0:0:0';
    my $test5 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=1	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:0:0:0:0	4:0:0:0:0:0:0:0:0:0';
    my $test6 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=1	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:0:0:0:0	0:4:0:0:0:0:0:0:0:0';
    my $test7 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=1	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:0:0:0:0	5:1:1:0:0:0:0:0:0:0';
     my $test8 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=1	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:0:0:0:0	0:5:0:0:0:0:0:0:0:0';
     my $test9 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=1	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:0:0:0:0	4:1:0:0:0:0:0:0:0:0';
    my $test10 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=1	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:0:0:0:0	1:4:0:0:0:0:0:0:0:0';
    my $test11 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=2	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:0:0:0:0	0:5:0:0:0:0:0:0:0:0';
    my $test12 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=2	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:0:0:0:0	4:1:0:0:0:0:0:0:0:0';
    my $test13 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=2	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:0:0:0:0	1:4:0:0:0:0:0:0:0:0';

    my $test14 = '22	16404839	.	G	GA	.	.	PC=I;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=0	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:0:0:0:0	1:4:0:0:0:0:0:0:0:0';
    my $test15 = '22	16404839	.	G	GA	.	.	PC=I;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=1	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:0:0:0:0	1:4:0:0:0:0:0:0:0:0';


    my $sub = $filter_hash->{test};

    my($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF);
    $CHROM  = [split("\t",$test1)]->[0];
    $POS    = [split("\t",$test1)]->[1];
    $FAIL   = 1;
    $PASS   = 0;
    #$VCF    = $$opts{vcf};

    is($filter_hash->{tag}, 'INFO/REP',"_test_F016 check the correct info tag has been set for the rule");

    $RECORD = [split("\t",$test1)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F016 pMt == 5");
    $RECORD = [split("\t",$test2)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F016 pMt == 5 with different pos/neg ratios");
    $RECORD = [split("\t",$test3)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F016 pMt == 5 with different pos/neg ratios");
    $RECORD = [split("\t",$test4)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F016 pMt == 6");
    $RECORD = [split("\t",$test5)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_F016 pMt == 4");
    $RECORD = [split("\t",$test6)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_F016 pMt == 4 with different pos/neg ratios");
    $RECORD = [split("\t",$test7)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F016 pMt == 5 bMt > 0");

    $MATCH = 1;
    $RECORD = [split("\t",$test8)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_F016 pMt == 5 bMt > 0 rep == 1");
    $MATCH = 1;
    $RECORD = [split("\t",$test9)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F016 pMt == 5 bMt > 0 rep == 1 pMtPos > 1 pMtPos == 0");
    $MATCH = 1;
    $RECORD = [split("\t",$test10)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F016 pMt == 5 bMt > 0 rep == 1 pMtPos == 0 pMtNeg > 1");

    $MATCH = 2;
    $RECORD = [split("\t",$test11)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_F016 pMt == 5 bMt > 0 rep == 2: pindel calls in only one direction");
    $MATCH = 2;
    $RECORD = [split("\t",$test12)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_F016 pMt == 5 bMt > 0 rep == 2: pindel calls in both directions");
    $MATCH = 2;
    $RECORD = [split("\t",$test13)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_F016 pMt == 5 bMt > 0 rep == 2: pindel calls in both direction");

    $MATCH = 0;
    $RECORD = [split("\t",$test14)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F016 pMt == 5 bMt > 0 rep == 0: insertion pindel calls in both directions");
    $MATCH = 1;
    $RECORD = [split("\t",$test15)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_F016 pMt == 5 bMt > 0 rep == 1: insertion pindel calls in both directions");
  };
}

sub _test_F017{
	#fail; #simple repeat filter
}

sub _test_F018{
	my ($filter_hash) = @_;
  subtest "Test rule F018" => sub {
    my $test1 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:5:5:0:0	0:0:0:0:0:0:5:5:0:0';
    my $test2 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:6:4:0:0	0:0:0:0:0:0:5:5:0:0';
    my $test3 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:5:5:0:0	0:0:0:0:0:0:6:4:0:0';
    my $test4 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:5:4:0:0	0:0:0:0:0:0:5:5:0:0';
    my $test5 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:5:5:0:0	0:0:0:0:0:0:5:4:0:0';
    my $test6 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:5:6:0:0	0:0:0:0:0:0:5:5:0:0';
    my $test7 = '22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	0:0:0:0:0:0:5:5:0:0	0:0:0:0:0:0:5:6:0:0';

    my $sub = $filter_hash->{test};

    my($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF);
    $CHROM  = [split("\t",$test1)]->[0];
    $POS    = [split("\t",$test1)]->[1];
    $FAIL   = 1;
    $PASS   = 0;
    #$VCF    = $$opts{vcf};

    is($filter_hash->{tag}, 'INFO/LEN',"_test_F018 check the correct info tag has been set for the rule");

    $RECORD = [split("\t",$test1)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F014 pWt == 10 pMt == 10");
    $RECORD = [split("\t",$test2)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F014 pWt == 10 pMt == 10 with different pos/neg ratios");
    $RECORD = [split("\t",$test3)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F014 pWt == 10 pMt == 10 with different pos/neg ratios");
    $RECORD = [split("\t",$test4)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_F014 pWt == 10 pMt == 9");
    $RECORD = [split("\t",$test5)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $FAIL,"_test_F014 pWt == 9 pMt == 10");
    $RECORD = [split("\t",$test6)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F014 pWt == 11 pMt == 10");
    $RECORD = [split("\t",$test7)];
    is($sub->($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF), $PASS,"_test_F014 pWt == 10 pMt == 11");
  };
}
