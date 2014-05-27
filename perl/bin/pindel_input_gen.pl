#!/usr/bin/perl

BEGIN {
  use Cwd qw(abs_path);
  use File::Basename;
  unshift (@INC,dirname(abs_path($0)).'/../lib');
};

use strict;
use warnings FATAL => 'all';
use autodie qw(:all);

use File::Path qw(remove_tree make_path);
use Getopt::Long;
use File::Spec;
use Pod::Usage qw(pod2usage);
use List::Util qw(first);
use Const::Fast qw(const);
use File::Temp;

use Sanger::CGP::Pindel::InputGen;

{
  my $options = setup();

  my $generator = Sanger::CGP::Pindel::InputGen->new($options->{'bam'});
  $generator->set_threads($options->{'threads'});
  $generator->set_outdir($options->{'outdir'});
  $generator->run;
}

sub setup {
  my %opts;
  GetOptions( 'h|help' => \$opts{'h'},
              'm|man' => \$opts{'m'},
              't|threads=i' => \$opts{'threads'},
              'o|outdir=s' => \$opts{'outdir'},
              'b|bam=s' => \$opts{'bam'},
  ) or pod2usage(2);

  pod2usage(-verbose => 2) if(defined $opts{'h'});
  pod2usage(-verbose => 1) if(defined $opts{'m'});

  # then check for no args:
  my $defined;
  for(keys %opts) { $defined++ if(defined $opts{$_}); }
  pod2usage(-msg  => "\nERROR: Options must be defined.\n", -verbose => 2,  -output => \*STDERR) unless($defined);

  pod2usage(-msg  => "\nERROR: 'outdir' must be defined.\n", -verbose => 2,  -output => \*STDERR) unless(defined $opts{'outdir'});
  pod2usage(-msg  => "\nERROR: 'bam' must be defined.\n", -verbose => 2,  -output => \*STDERR) unless(defined $opts{'bam'});

  $opts{'threads'} ||= 1;

  die "'outdir' is not writable: $opts{outdir}\n" if(-e $opts{'outdir'} && !-w _);
  unless(-e $opts{'outdir'}) {
    make_path $opts{'outdir'} or die "Failed to create $opts{outdir}";
  }

  return \%opts;
}

__END__

=head1 pindel_input_gen.pl

Generate input files for pindel from a BAM file.

=head1 SYNOPSIS

pindel_input_gen.pl [options]

  Required parameters:
    -bam       -b   BAM file to process.
    -outdir    -o   Folder to output result to.
    -threads   -t   Number of threads to use. [1]

  Other:
    -help      -h   Brief help message.
    -man       -m   Full documentation.

  File list can be full file names or wildcard, e.g.
    pindel_input_gen.pl -t 16 -o myout -b some.bam

  MEMORY:   Requires ~1GB of RAM.
            Linear growth with threads.

  RUNTIME:  ~100 sec. per million pairs.
            ~50 sec. when 2 threads used.
            Greater than 2 threads not recommended

=head1 OPTIONS

=over 8

=item B<-bam>

Mapped BAM file to be prepared.

=item B<-outdir>

Output location, one file will be generated per ref-seq with anchors.

Temporary files will also be written to this area.

=item B<-threads>

Number of threads to be used in processing.  2 is recommended, but memory usage is ~doubled.

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=back
