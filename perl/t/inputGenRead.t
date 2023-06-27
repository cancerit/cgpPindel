

use strict;
use Test::More;
use Test::Fatal;
use Const::Fast qw(const);

const my $MODULE => 'Sanger::CGP::Pindel::InputGen::Read';

const my $READ => <<'READ';
HS16_08993:1:2314:21033:57219#78	147	1	1223958	12	10M2I63M	=	1223723	-308	TCTGCCCCGTCCCCCGTGTCTCTGCCCCGTCCCCCGTGTCTCTGCTCCGTCCTCCGTGTCTCTGCCCCGTCCCGT	*,+++'A>>(B6F5B>54CB5+,?BE6@47?HDDEE>AE,AH?F7FE@F@7B,DADFHGEEGGHFGF@FFGCDD?	MD:Z:5T17T26C22	RG:Z:1141553	XG:i:2	AM:i:12	NM:i:5	SM:i:12	XM:i:3	XO:i:1	XT:A:M
READ

const my $EXP_FRAC_PBQ => 0.173333333333333;

my $obj;
subtest 'Initialisation checks' => sub {
  use_ok($MODULE);
};

my $cleaned = $READ;
chomp $cleaned;
my ($rg) = $cleaned =~ m/\tRG:Z:([^\t]+)/;
$obj = new_ok($MODULE, [\$cleaned, 2]);

is($obj->frac_pbq_poor, $EXP_FRAC_PBQ, 'Check poor qual PBQ fraction');

done_testing();
