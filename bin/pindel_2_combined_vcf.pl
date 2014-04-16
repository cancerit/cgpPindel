#! /usr/bin/perl

BEGIN {
  use Cwd qw(abs_path);
  use File::Basename;
  push (@INC,dirname(abs_path($0)).'/../lib');
};

use strict;
use warnings FATAL => 'all';
use autodie qw(:all);

use Getopt::Long;
use Try::Tiny;
use Carp;
use Pod::Usage qw(pod2usage);
use Data::Dumper;

use Bio::DB::Sam;

use Sanger::CGP::Pindel::OutputGen::CombinedRecordGenerator;
use Sanger::CGP::Pindel::OutputGen::SequentialIdGenerator;
use Sanger::CGP::Pindel::OutputGen::UuIdGenerator;
use Sanger::CGP::Pindel::OutputGen::VcfConverter;

use Sanger::CGP::Vcf::Contig;
use Sanger::CGP::Vcf::Sample;
use Sanger::CGP::Vcf::BamUtil;

{
  my $opts = setup();
  
  my $output;
  my $output_location = '';
  my $current_path = '';
  my @file_paths = ();
  
  if($opts->{input}){
    push(@file_paths, $opts->{input});
  }
  	
  if($opts->{inputdir}){
    opendir(my $dh, $opts->{inputdir}) or die "Cannot open |".$opts->{d}."| for reading: $!";
    push(@file_paths, map {join('/',$opts->{inputdir},$_)} grep {m/.*_((D)|(SI))$/} readdir $dh);
    closedir $dh;
  }
  	
  my $wt_sam = Bio::DB::Sam->new(-bam => $opts->{wt}, -fasta => $opts->{r});
  my $mt_sam = Bio::DB::Sam->new(-bam => $opts->{mt}, -fasta => $opts->{r});

  ## combine the header text to get a combined list of contigs.
  my $contigs    = Sanger::CGP::Pindel::OutputGen::BamUtil::parse_contigs($wt_sam->header->text . $mt_sam->header->text);
  my $wt_samples = Sanger::CGP::Pindel::OutputGen::BamUtil::parse_samples($wt_sam->header->text,$opts->{wts},$opts->{wtp},$opts->{wta},$opts->{a},'Wild type');
  my $mt_samples = Sanger::CGP::Pindel::OutputGen::BamUtil::parse_samples($mt_sam->header->text,$opts->{mts},$opts->{mtp},$opts->{mta},$opts->{a},'Mutant');

  ## make sure we have the correct number of samples in the bam files.
  die "No samples found in wild type bam file." if(scalar values %$wt_samples == 0);
  die "No samples found in mutant type bam file." if(scalar values %$mt_samples == 0);
  die "Multiple samples found in wild type bam file." if(scalar values %$wt_samples > 1);
  die "Multiple samples found in mutant type bam file." if(scalar values %$mt_samples > 1);
 
  my $wt_sample = (values(%$wt_samples))[0];
  my $mt_sample = (values(%$mt_samples))[0];
  my $fai = Bio::DB::Sam::Fai->load($opts->{r});

  my $wt_out_sam_path;
  my $mt_out_sam_path;
  my $wt_out_sam_fh;
  my $mt_out_sam_fh;

  if(defined $opts->{'samoutput'}){
  	my $base_path = $opts->{'samoutput'} =~ s/\..+$//r;
    $wt_out_sam_path = $base_path . '_wt.sam';
    $mt_out_sam_path = $base_path . '_mt.sam';

    open($wt_out_sam_fh, ">", $wt_out_sam_path) or die "Cannot open |$wt_out_sam_path| for writing: $!";
    open($mt_out_sam_fh, ">", $mt_out_sam_path) or die "Cannot open |$mt_out_sam_path| for writing: $!";
  }

  my $record_converter = new Sanger::CGP::Pindel::OutputGen::VcfConverter(
    -contigs => [values %$contigs]
  );

  my $id_gen;
  if(defined $opts->{'g'}){
    $id_gen = new Sanger::CGP::Pindel::OutputGen::SequentialIdGenerator(-start => $opts->{'g'});
  }else{
  	$id_gen = new Sanger::CGP::Pindel::OutputGen::UuIdGenerator();
  }

  try{
    if($opts->{'output'}){
      open($output, '>', $opts->{'output'}) or die "Cannot open |".$opts->{'output'}."| for reading: $!";
      $output_location = $opts->{'output'};
    }else{
      $output = \*STDOUT;
      $output_location = 'STDOUT';
    }

    ## loop though all the pindel files for the given locations.
    if(scalar @file_paths){
      my $print_header = 1;
      foreach $current_path (@file_paths){
      	open(my $fh, "<", $current_path) or croak "Cannot open |$current_path| for reading: $!";
        _process_fh($fh, $output, $print_header, $wt_out_sam_fh, $mt_out_sam_fh, $opts, $record_converter, $fai, $wt_sam, $mt_sam, $wt_sample, $mt_sample, $id_gen, $opts->{r});
        close($fh) or croak "Unable to close |$current_path| for reading: $!";
        $print_header = 0;
      }
    }else{
    	$current_path = 'STDIN';
    	_process_fh(\*STDIN, $output, 1, $wt_out_sam_fh, $mt_out_sam_fh, $opts, $record_converter, $fai, $wt_sam, $mt_sam, $wt_sample, $mt_sample, $id_gen, $opts->{r});
    }

  }catch{
    die 'Error! Reading: |'.$current_path."| Writing to: |$output_location|\n$_";
  }finally{
    close $output or die "Unable to close |".$opts->{'output'}."|: $!" if($opts->{'output'});
    close $wt_out_sam_fh or die "Unable to close $wt_out_sam_path|: $!" if(defined $wt_out_sam_fh);
    close $mt_out_sam_fh or die "Unable to close $mt_out_sam_path|: $!" if(defined $mt_out_sam_fh);
  }
}

sub _process_fh{
  my($fh, $output, $header, $wt_out_sam_fh, $mt_out_sam_fh, $opts, $record_converter, $fai, $wt_sam, $mt_sam, $wt_sample, $mt_sample, $id_gen, $ref) = @_;

  my $record_generator = new Sanger::CGP::Pindel::OutputGen::CombinedRecordGenerator(
    -wt_sam => $wt_sam,
    -mt_sam => $mt_sam,
    -mutant_sample_name => $mt_sample->name,
    -fai => $fai,
    -fh => $fh
  );
  

  if($header){
    ## write the header... we need to test if S2 is preant first...
    my $source = basename($0). '_v'. Sanger::CGP::Pindel::OutputGen->VERSION;
    my $include_s2 = $record_generator->version eq 'v01' ? 1 : 0;
    
    my @process_logs = ();
    
    push @process_logs, new Sanger::CGP::Vcf::VcfProcessLog(
	  -input_vcf_source => 'Pindel',
      -input_vcf_ver => $record_generator->version,
	);
	
	push @process_logs, new Sanger::CGP::Vcf::VcfProcessLog(
	  -input_vcf_source => basename($0),
      -input_vcf_ver => Sanger::CGP::Pindel::OutputGen->VERSION,
      -input_vcf_param => $opts,
	);
    
    print $output $record_converter->gen_header($wt_sample,$mt_sample, \@process_logs, $include_s2, $ref, $source) or croak("Unable to write VCF header: $!");
  }
  
  my ($active_sam_fh, $sample, $strand);
  while(my $record = $record_generator->next_record){

    $record->id($id_gen->next);

    #unless($opts->{s} && ($record->p_wt_pos + $record->p_wt_neg) > ($record->p_mt_pos + $record->p_mt_neg)){
    #unless($opts->{s} && $record->valid == 0){
   	
    unless(
            $opts->{s} &&
            (
              $record->valid == 0 ||
              (
                ($record->p_mt_pos + $record->p_mt_neg) < 2 #&&
                #($record->p_wt_pos + $record->p_wt_neg) < 2
              )
            )
          ){
    	
      ## write the record
      print $output $record_converter->gen_record($record) or croak("Unable to write VCF record: $!") ;

      ## only write to sam files if they have been defined.
      ## these can be used for viewing in g/jbrowse.
      if($mt_out_sam_fh){
        foreach $sample (keys %{$record->reads}){

          $active_sam_fh = $sample eq $mt_sample->name ? $mt_out_sam_fh : $wt_out_sam_fh;

          for $strand ('+','-'){
            foreach my $read_arr (@{$record->reads->{$sample}->{$strand}}){
              print $active_sam_fh join("\t",@$read_arr)."\n" or die "Unable to write sam line: $!";
            }
          }

        }
      }
    }
  }
}

sub setup{
  my %opts;
  my @random_args;
  GetOptions( 'h|help' => \$opts{'h'},
              'm|man' => \$opts{'m'},
              'v|version' => \$opts{'v'},
              'i|input=s' => \$opts{'input'},
              'd|inputdir=s' => \$opts{'inputdir'},
              'so|samoutput=s' => \$opts{'samoutput'},
              'o|output=s' => \$opts{'output'},
              'r|ref=s' => \$opts{'r'},
              'w|wtbam=s' => \$opts{'wt'},
              'm|mtbam=s' => \$opts{'mt'},
              'ws|wtstudy=s' => \$opts{'wts'},
              'ms|mtstudy=s' => \$opts{'mts'},
              'wp|wtprot=s' => \$opts{'wtp'},
              'mp|mtprot=s' => \$opts{'mtp'},
              'wa|wtacc=s' => \$opts{'wta'},
              'ma|mtacc=s' => \$opts{'mta'},
              'a|accsource=s' => \$opts{'a'},
              's|skipwt' => \$opts{'s'},
              'g|idstart=i' => \$opts{'g'},
              '<>' => sub{push(@random_args,shift(@_));}
  ) or pod2usage(2);

  my $version = Sanger::CGP::Pindel::OutputGen->VERSION;

  if(defined $opts{'v'}){
    print "Version: $version\n";
    exit;
  }

  pod2usage(-verbose => 1) if(defined $opts{'h'});
  pod2usage(-verbose => 2) if(defined $opts{'m'});

  pod2usage(-message  => "\nERROR: unrecognised commandline arguments: ".join(', ',@random_args).".\n", -verbose => 1,  -output => \*STDERR) if(scalar @random_args) ;

  pod2usage(-message  => "\nERROR: w|wtbam must be defined.\n", -verbose => 1,  -output => \*STDERR) unless($opts{'wt'});
  pod2usage(-message  => "\nERROR: w|wtbam |".$opts{'wt'}."| must be a valid file.\n", -verbose => 1,  -output => \*STDERR) unless(-f $opts{'wt'});
  pod2usage(-message  => "\nERROR: w|wtbam |".$opts{'wt'}."| cannot locate index file.\n", -verbose => 1,  -output => \*STDERR) unless(-f $opts{'wt'}.'.bai');

  pod2usage(-message  => "\nERROR: m|mtbam must be defined.\n", -verbose => 1,  -output => \*STDERR) unless($opts{'mt'});
  pod2usage(-message  => "\nERROR: m|mtbam |".$opts{'mt'}."| must be a valid file.\n", -verbose => 1,  -output => \*STDERR) unless(-f $opts{'mt'});
  pod2usage(-message  => "\nERROR: w|wtbam |".$opts{'mt'}."| cannot locate index file.\n", -verbose => 1,  -output => \*STDERR) unless(-f $opts{'mt'}.'.bai');

  pod2usage(-message  => "\nERROR: r|ref must be defined.\n", -verbose => 1,  -output => \*STDERR) unless($opts{'r'});
  pod2usage(-message  => "\nERROR: r|ref |".$opts{'r'}."| must be a valid file.\n", -verbose => 1,  -output => \*STDERR) unless(-f $opts{'r'});

  return \%opts;
}


__END__

=head1 NAME

Pindel2CombinedVcf.pl - Takes a raw Pindel file and a pair of wild type/mutant bam files and produces a combined .vcf file.

=head1 SYNOPSIS

Pindel2CombinedVcf.pl [options]

  Required parameters:
    -output    -o   File path to output. Defaults to STDOUT.
    -samoutput -so  File path-stub to sam file output. If presant will create two sam files <path-stub>_wt|mt.sam containing the pindel call reads.
    -wtbam     -wt  File path to the wild type bam file (index must be at the same location).
    -mtbam     -wt  File path to the mutant bam file (index must be at the same location).
    -ref       -r   File path to the reference file used to provide the coordinate system.
    
  Optional parameters:
    
    -input     -i   File path to read in.
    -inputdir  -d   Directory path to read in. All files ending in _D or _SI.
    
    (If neither -i nor -d are set, input will be read from STDIN)
    
    -wtstudy   -wts String representing the wild type sample study.
    -mtstudy   -mts String representing the mutant sample study.
    -wtprot    -wts String representing the wild type sample sequence protocol (e.g. genomic, targeted, RNA-seq).
    -mtprot    -mts String representing the mutant sample sequence protocol (e.g. genomic, targeted, RNA-seq).
    -wtacc     -wta String representing the wild type sample accession id.
    -mtacc     -mta String representing the mutant sample accession id.
    -accsource -a   String representing the source of the accession id (e.g. COSMIC_SAMPLE_ID, EGA etc...)
    -skipwt    -s   If presant, will skip variants where there are more wt calls than mt.
    -idstart   -g   Will set a sequential id generator to the given integer value. If not presant will assign UUIDs to each variant.

  Other:
    -help      -h   Brief help message.
    -man       -m   Full documentation.
    -version   -v   Prints the version number.

    Pindel2CombinedVcf.pl -i my.bam -o my.bam.bas -wt wt.bam -mt mt.bam

=head1 OPTIONS

=over 8

=item B<-input>

File path to read. Accepts only raw pindel files.

=item B<-samoutput>

File path-stub to sam file output. If presant will create two sam files <path-stub>_wt|mt.sam containing the pindel call reads. These files are used to display the pindel reads in various browsers.

=item B<-output>

File path to output data. If this option is omitted the script will attempt to write to STDOUT. The

=item B<-wtbam>

File path to read. Accepts only .bam files. This file is treated as the wild type sample.

=item B<-mtbam>

File path to read. Accepts only .bam files. This file is treated as the mutant sample.

=item B<-wtstudy>

String identifying the study to which the wild type sample belongs.

=item B<-mtstudy>

String identifying the study to which the mutant sample belongs.

=item B<-wtprot>

String identifying the sequence protocol that wild type sample was sequenced under (e.g. genomic, targeted, RNA-seq).

=item B<-mtprot>

String identifying the sequence protocol that mutant sample was sequenced under (e.g. genomic, targeted, RNA-seq).

=item B<-wtacc>

String identifying the wild type sample accession identifier.

=item B<-mtacc>

String identifying the mutant sample accession identifier.

=item B<-accsource>

String identifying the accession source (e.g. COSMIC_SAMPLE_ID, EGA etc...).

=item B<-skipwt>

If presant, will skip variants where there are more wt calls than mt. Unfortunatelty the variant will still be processed in the background.

=item B<-idstart>

Will set the id generator to the given integer value. If not presant will assign UUIDs to each variant.

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=item B<-version>

Prints the version number and exits.

=back

=head1 DESCRIPTION

B<Pindel2CombinedVcf.pl> will attempt to generate a vcf file from a Pindel output file and bwa pileup information.
For evey variant called by Pindel a pileup will be performed and the results merged into a single vcf record.

=cut