use strict;
use File::Temp qw(tempdir);
use File::Path qw(make_path);
use Test::More;
use Test::Fatal;
use File::Spec::Functions;
use Const::Fast qw(const);
use FindBin qw($Bin);

const my $MODULE => 'Sanger::CGP::Pindel::OutputGen::VcfBlatAugment';
const my $DATA => "$Bin/data/blat";
const my @HEADER_ENDS => do {
    no warnings 'qw';
    qw(#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT PD26988a);
};
const my $HEADER_LINES => 3393;

const my $DATA_ARR_D => [qw(chr10 11201 id CA C 390 . PC=D;RS=11201;RE=11205;LEN=1;S1=11;S2=849.236;REP=3 GT:PP:NP ./.:10:0)];
const my $RES_ARR_D => [qw(chr10 11201 id CA C 390 . PC=D;RS=11201;RE=11205;LEN=1;S1=11;S2=849.236;REP=3 GT:PP:NP ./.:10:0:12:3:0.015:8:3:0.006:0.423)];
const my $DATA_ARR_DI => [qw(chr10 22777 id AGAAACTGTG ACTGTGAGATAGATATATATAGATAGATATAT 105 . PC=DI;RS=22777;RE=22787;LEN=9;S1=6;REP=0 GT:PP:NP ./.:0:5)];
const my $RES_ARR_DI => [qw(chr10 22777 id AGAAACTGTG ACTGTGAGATAGATATATATAGATAGATATAT 105 . PC=DI;RS=22777;RE=22787;LEN=9;S1=6;REP=0 GT:PP:NP ./.:0:5:2:2:0.002:2:8:0.017:0.714)];
const my $DATA_ARR_SI => [qw(chr10 11643 id C CG 150 . PC=I;RS=11643;RE=11649;LEN=1;S1=6;S2=421.908;REP=4 GT:PP:NP ./.:0:5)];
const my $RES_ARR_SI => [qw(chr10 11643 id C CG 150 . PC=I;RS=11643;RE=11649;LEN=1;S1=6;S2=421.908;REP=4 GT:PP:NP ./.:0:5:7:10:0.014:7:10:0.005:0.500)];

my ($stdout_fh, $buffer);
my ($sam_stdout_fh, $sam_buffer);

subtest 'Initialisation checks' => sub {
  use_ok($MODULE);
  new_vba(catdir($DATA, 'D'));
};

subtest 'Header checks' => sub {
  my $vba = new_vba(catdir($DATA, 'D'));
  ok($vba->output_header);
  my @lines = split /\n/, $buffer;
  is(scalar @lines, $HEADER_LINES, 'Expected number of header lines');
  is($lines[-1], join("\t", @HEADER_ENDS), 'Expected final header line');
};

subtest 'Simple Deletion checks' => sub {
  my $vba = new_vba(catdir($DATA, 'D'));
  my @tmp = @{$DATA_ARR_D};
  $vba->blat_record(\@tmp);
  is_deeply(\@tmp, $RES_ARR_D);
};

subtest 'Simple Insertion checks' => sub {
  my $vba = new_vba(catdir($DATA, 'SI'));
  my @tmp = @{$DATA_ARR_SI};
  $vba->blat_record(\@tmp);
  is_deeply(\@tmp, $RES_ARR_SI);
};

subtest 'Complex event checks' => sub {
  my $vba = new_vba(catdir($DATA, 'DI'));
  my @tmp = @{$DATA_ARR_DI};
  $vba->blat_record(\@tmp);
  print join q{ }, @tmp;
  is_deeply(\@tmp, $RES_ARR_DI);
};

done_testing();


sub new_vba {
  my $dir = shift;
  my $tmp = '/tmp/pindel_test_stuff';
  make_path($tmp);
  my $obj = new_ok($MODULE, [
    input => catfile($dir, 'test.vcf'),
    ref => catfile($DATA, 'chr10_1-23700.fa'),
    ofh => buffer_fh(),
    outpath => $tmp,
    hts_files => [catfile($DATA, 'test.bam')],
  ]);
  my $sample = $obj->{vcf_sample_order}->[0];
  $obj->{sfh}->{$sample} = sam_buffer_fh();
  unlink $tmp;
  return $obj;
}

sub buffer_fh {
  if(defined $stdout_fh) {
    close $stdout_fh;
  }
  $buffer = q{};
  open $stdout_fh, ">", \$buffer or die $!;
  return $stdout_fh;
}

sub sam_buffer_fh {
  if(defined $sam_stdout_fh) {
    close $sam_stdout_fh;
  }
  $sam_buffer = q{};
  open $sam_stdout_fh, ">", \$sam_buffer or die $!;
  return $sam_stdout_fh;
}



__END__
subtest 'corrupt_pindel_input checks' => sub {
  # File of correct size (no NUL)
  is( Sanger::CGP::Pindel::InputGen::corrupt_pindel_input("$DATA/inputGen-goodfile.txt.gz", 22),
      undef,
      'corrupt_pindel_input - well formed compressed file, expected size');
  # File of incorrect size (no NUL)
  is( Sanger::CGP::Pindel::InputGen::corrupt_pindel_input("$DATA/inputGen-goodfile.txt.gz", 9),
      "$DATA/inputGen-goodfile.txt.gz",
      'corrupt_pindel_input - well formed compressed file, UNexpected size');
  # File of incorrect size + NUL
  is( Sanger::CGP::Pindel::InputGen::corrupt_pindel_input("$DATA/inputGen-NUL.txt.gz", 22),
      "$DATA/inputGen-NUL.txt.gz",
      'corrupt_pindel_input - NUL character in compressed file, expected size');
  # File of correct size + NUL
  is( Sanger::CGP::Pindel::InputGen::corrupt_pindel_input("$DATA/inputGen-NUL.txt.gz", 9),
      "$DATA/inputGen-NUL.txt.gz",
      'corrupt_pindel_input - NUL character in compressed file, UNexpected size');
};

subtest 'reads_to_disk checks' => sub{
  $obj = new_ok($MODULE,
                [ "$DATA/test.bam",
                  undef,
                  "$DATA/genome_22.fa"]);
  my $out_folder = tempdir( 'pindelTests_XXXX', CLEANUP => 1 );

  $obj->set_outdir($out_folder);
  ok($obj->reads_to_disk($RECORD_SET), 'create output');
  is($obj->{'rname_bytes'}->{'22'}, $RECORD_OUT_BYTES, 'Verify bytes written captured');

  ok($obj->validate, 'validate returns true')
};

done_testing();
