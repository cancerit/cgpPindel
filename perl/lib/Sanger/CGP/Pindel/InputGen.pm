package Sanger::CGP::Pindel::InputGen;

########## LICENCE ##########
# Copyright (c) 2014 Genome Research Ltd. 
#  
# Author: Keiran Raine <cgpit@sanger.ac.uk> 
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
use English qw( -no_match_vars );
use autodie qw(:all);
use warnings FATAL => 'all';
use Config; # so we can see if threads are enabled
use Const::Fast qw(const);
use Carp qw( croak );
use File::Which qw(which);
use Try::Tiny qw(try catch finally);
use File::Spec;
use Data::Dumper;

use Sanger::CGP::Pindel;

# add signal handler for any interrupt so that cleanup of temp is handled correctly.
use File::Temp;
$SIG{INT} = sub { exit };
$SIG{TERM} = sub { exit };

use threads;

use Sanger::CGP::Pindel::InputGen::SamHeader;
use Sanger::CGP::Pindel::InputGen::PairToPindel;
use Sanger::CGP::Pindel::InputGen::Pair;

const my $PAIRS_PER_THREAD => 500_000;

const my $BAMCOLLATE => q{%s outputformat=sam colsbs=268435456 collate=1 classes=F,F2 exclude=DUP,SECONDARY,QCFAIL,SUPPLEMENTARY T=%s filename=%s};

sub new {
  my ($class, $bam) = @_;
  my $self = {'rname_fhs' => {},
              'threads' => 1, };
  bless $self, $class;
  $self->set_input($bam) if(defined $bam);
  return $self;
}

sub set_input {
  my ($self, $bam) = @_;
  croak "set_input requires a value for 'bam'" unless(defined $bam);
  die "Does not appear to be a BAM file: $bam" unless($bam =~ m/\.bam$/);
  die "File does not exist : $bam" unless(-e $bam);
  die "File appears to be empty : $bam" unless(-s _);
  $self->{'bam'} = $bam;
  return $self->{'bam'};
}

sub set_threads {
  my ($self, $threads) = @_;
  croak "set_threads requires a value for 'threads'" unless(defined $threads);
  die "Number of threads should be a positive integer: $threads" if($threads !~ m/^[[:digit:]]+$/ || $threads < 1);

  $self->{'threads'} = $threads;
  return $self->{'threads'};
}

sub set_outdir {
  my ($self, $outdir) = @_;
  croak "set_output requires a value for 'outdir'" unless(defined $outdir);
  die "Output folder does not exist: $outdir" unless(-e $outdir);
  die "Output folder is not writeable: $outdir" unless(-w _);
  $self->{'outdir'} = $outdir;
  return $self->{'outdir'};
}

sub run {
  my $self = shift;

  # set up values used in loop;
  my $pair_count = 0;
  my $threads = $self->{'threads'};

  # setup the temp dir:
  my $tmpdir = File::Temp->newdir( File::Spec->catdir($self->{'outdir'}, 'tmpXXXX') );

  my $collate = which('bamcollate2');

  # ensure that commands containing pipes give appropriate errors
  my $command .= sprintf $BAMCOLLATE, $collate, File::Spec->catfile($tmpdir, 'collate_tmp'), $self->{'bam'};
    try {
    my ($rg_pis, $sample_name);
    my $head_ob = Sanger::CGP::Pindel::InputGen::SamHeader->new($self->{'bam'});
    my $collate_start = time;
    my ($pid, $process);
    $pid = open $process, q{-|}, $command or croak 'Could not fork: '.$OS_ERROR;

    my @this_batch;
    while (my $tmp = <$process>) {
      if(defined $rg_pis) {
        ## SAM processing
        # get second in pair
        my $tmp_b = <$process>;
        chomp $tmp;
        chomp $tmp_b;
        push @this_batch, $tmp, $tmp_b;
        $pair_count++;
        if($pair_count == $PAIRS_PER_THREAD) {
          my $collate_took = time - $collate_start;
          warn "Collated $pair_count readpairs (in $collate_took sec.)\n";
          $self->_process_set($rg_pis, $sample_name, \@this_batch);
          @this_batch = ();
          $pair_count = 0;
          $collate_start = time;
        }

        next;
      }
      if($tmp =~ m/^\@/) {
        $head_ob->add_header_line($tmp);
      }
      else {
        ($rg_pis, $sample_name) = $head_ob->rgs_to_offsets_and_sample_name;
        redo; # re-process this line of input now header object is defined
      }
    }
    close $process;
    my $collate_took = time - $collate_start;
    warn "Collated $pair_count readpairs (${collate_took}s)\n";

    $self->_process_set($rg_pis, $sample_name, \@this_batch);
    $self->_completed_threads if($threads > 1);
  }
  catch {
    die qq{An error occurred while running:\n\t$command\nERROR: $_};
  };
}

sub _process_set {
  my ($self, $rg_pis, $sample_name, $pairs) = @_;
  my $max_threads = $self->{'threads'};
  if($max_threads > 1) {
    my $existing_threads = threads->list(threads::all);
    if(@{$pairs} == 0) {
      $self->_completed_threads;
      return 1;
    }
    # first clear all old processing
    if($existing_threads == $max_threads) {
      sleep 1 while(threads->list(threads::joinable) < $max_threads);
      $self->_completed_threads;
    }
    # start new thread
    my ($thr) = threads->create(\&reads_to_pindel, $existing_threads, $rg_pis, $sample_name, @{$pairs});
  }
  else {
    $self->reads_to_disk([reads_to_pindel(-1,  $rg_pis, $sample_name, $pairs)]);
  }
}

sub _completed_threads {
  my $self = shift;
  my @output_objects;
  my $all_threads = threads->list(threads::all);
  sleep 1 while(threads->list(threads::joinable) != $all_threads);
  for my $thr(threads->list(threads::joinable)) {
    my @data = $thr->join;
    if(my $err = $thr->error) { die "Converter thread error: $err\n"; }
    push @output_objects, \@data;
  }
  $self->reads_to_disk(\@output_objects);
}

sub reads_to_pindel {
  my ($thread, $rg_pis, $sample_name, @reads) = @_;
  my $tid = 0;
  $tid = threads->tid unless($thread == -1);
  warn "Thread Worker $tid: started\n";

  my $convert_start = time;
  @reads = @{$reads[0]} if(ref $reads[0] eq 'ARRAY');

  my $total_reads = scalar @reads;
  my @records;
  my $to_pindel = Sanger::CGP::Pindel::InputGen::PairToPindel->new($sample_name, $rg_pis);
  my $total_pairs = $total_reads / 2;
  for(1..$total_pairs) {
    my $pair = Sanger::CGP::Pindel::InputGen::Pair->new(\shift @reads, \shift @reads);
    next unless($pair->keep_pair);
    $to_pindel->set_pair($pair);
    push @records, @{$to_pindel->pair_to_pindel};

  }
  my $retained = scalar @records;
  my $excluded = $total_pairs - $retained;
  my $convert_took = time - $convert_start;
  warn "Thread Worker $tid: Excluded $excluded/$total_pairs (${convert_took}s)\n";
  warn "Thread Worker $tid: Generated $retained records\n";
  return \@records if($thread == -1);
  return @records;
}

sub reads_to_disk {
  my ($self, $record_sets) = @_;
  my $start = time;
  my %grouped;
  for my $set (@{$record_sets}) {
    for my $rec (@{$set}) {
      my ($rname) = $rec =~ m/\n[\+-]\t([^\t]+)/;
      push @{$grouped{$rname}}, $rec;
    }
  }
  my $rname_fh = $self->{'rname_fhs'};
  for my $rname(keys %grouped) {
    my $mode = '>>';
    unless(exists $rname_fh->{$rname}) {
      my $path = File::Spec->catfile($self->{'outdir'}, $rname.'.txt');
      $rname_fh->{$rname} = $path;
      $mode = '>';
    }
    open my $fh, $mode, $rname_fh->{$rname};
    for my $record(@{$grouped{$rname}}) {
      print $fh (join "\n", $record),"\n";
    }
    close $fh;
  }
}

1;

__END__

=head1 Sanger::CGP::Pindel::InputGen

Generate pindel input from BAM.

=head2 Constructor

=over 4

=item Sanger::CGP::Pindel::InputGen->new()

Create an input file generator, 2 options available

  my $in_gen = Sanger::CGP::Pindel::InputGen->new();
  $in_gen->set_input($bam_file);

or

  my $in_gen = Sanger::CGP::Pindel::InputGen->new($bam);
  $in_gen->set_threads($thread_count);
  $in_gen->set_outdir($out_folder);
  $in_gen->run();

=back

=head2 Methods

=over 4

=item set_input

Expects valid BAM file as input, will check presence and extension only
and croak if inappropriate.  The magic number for BAM is the same as any
gzip file (\037\213) so not used.

  $in_gen->set_input($bam_file);

=item set_threads(int)

Set the maximum number of threads to be used in processing.
When set to 1 no threading will be used.

  $in_gen->set_threads($max_threads);

Default 1.

=item set_outdir

Set the output folder. An input file is generated for each ref-seq found in the input BAM.
The calling program has the responsibility of checking/creating the folder.

  $in_gen->set_outdir($out_folder);

Function will check presence and write status and croak if missing/not-writable.

=item run

Start processing the BAM file and generating input files for pindel.

=back
