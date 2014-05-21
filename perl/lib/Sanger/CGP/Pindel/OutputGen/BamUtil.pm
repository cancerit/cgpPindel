package Sanger::CGP::Pindel::OutputGen::BamUtil;

use Sanger::CGP::Pindel;
use Sanger::CGP::Vcf::Contig;
use Sanger::CGP::Vcf::Sample;

use strict;
use autodie qw(:all);
use Carp;
use Data::Dumper;
use Capture::Tiny qw(capture_stdout);
use File::Which qw(which);

use Const::Fast qw(const);

const my $PG_TEMPLATE => "\@PG\tID:%s\tPN:%s\tCL:%s\tPP:%s\tDS:%s\tVN:%s";
const my $BAMSORT => ' inputformat=%s index=1 md5=1 O=%s indexfilename=%s md5filename=%s I=';

sub sam_to_sorted_bam {
  my ($path_prefix, $sam_files) = @_;
  my $bam_file = $path_prefix.'.bam';
  my $bai_file = $bam_file.'.bai';
  my $md5_file = $bam_file.'.md5';
  my $command = which('bamsort');
  $command .= sprintf $BAMSORT, 'sam', $bam_file, $bai_file, $md5_file;
  $command .= join q{ I=}, @{$sam_files};
  system($command);
  return 1;
}

sub pindel_header {
  my ($bam, $samp_type, $parent_pg, $cmd, $man_species, $man_assembly) = @_;
  my $command = which('samtools');
  $command .= ' view -H '.$bam;
  my $header = capture_stdout { system($command); };
  my @h_lines = split /\n/, $header;
  my $last_pg;
  for my $idx(0..(scalar @h_lines)-1) {
    my $h_line = $h_lines[$idx];
    if($h_line =~ m/^\@SQ/) {
      $h_lines[$idx] .= "\tSP:$man_species" if(defined $man_species && $h_line !~ m/\tSP:/);
      $h_lines[$idx] .= "\tAS:$man_assembly" if(defined $man_assembly && $h_line !~ m/\tAS:/);
    }
    $last_pg = $idx if($h_line =~ m/^\@PG/);
  }

  if(defined $last_pg) {
    my ($last_id) = $h_lines[$last_pg] =~ m/\tID:([^\t]+)/;
    my $preceeding_id;
    if($parent_pg) {
      ($preceeding_id) = $parent_pg =~ m|\\tID:([^\\]+)\\t|;
      $parent_pg =~ s/\\t/\t/g;
      $parent_pg =~ s/\tPP:\./\tPP:$last_id/;
      push @h_lines, $parent_pg;
      ($last_id) = $h_lines[-1] =~ m/\tID:[^\t]+/;
    }
    push @h_lines, pg_from_caller('pindel_to_combined_vcf', 'Converts text output of pindel to VCF and BAM alinments', $VERSION, $cmd, $preceeding_id);
  }
  return join "\n", @h_lines;
}

sub pg_from_caller {
  my ($id, $desc, $version, $cmd, $pp_id) = @_;
  $cmd =~ s/\t/\\t/g;

  return sprintf $PG_TEMPLATE,  $id,
                                $0, # program name
                                $cmd,
                                $pp_id || q{.},
                                $desc,
                                $version;
}

sub parse_contigs{
  my($header_txt, $man_species, $man_assembly) = @_;

  my $contigs = {};

  foreach my $line (split(/\n/,$header_txt)){
    if($line =~ /^\@SQ/){
      my ($name) = $line =~ /SN:([^\t]+)/;
      my ($length) = $line =~ /LN:([^\t]+)/;
      my ($assembly) = $line =~ /AS:([^\t]+)/;
      my ($species) = $line =~ /SP:([^\t]+)/;

      if(defined $man_species) {
        if(defined $species) {
          die "ERROR: Manually entered species ($man_species) doesn't match BAM file header ($species)\n"
            if($man_species ne $species);
        }
        else { $species = $man_species; }
      }
      die "ERROR: Species must be defined, check options/BAM header\n" unless(defined $species);

      if(defined $man_assembly) {
        if(defined $assembly) {
          die "ERROR: Manually entered assembly ($man_assembly) doesn't match BAM file header ($assembly)\n"
            if($man_assembly ne $assembly);
        }
        else { $assembly = $man_assembly; }
      }
      die "ERROR: Assembly must be defined, check options/BAM header\n" unless(defined $assembly);

      my $contig = new Sanger::CGP::Vcf::Contig(
        -name => $name,
        -length => $length,
        -assembly => $assembly,
        -species => $species
      );

      if(exists $contigs->{$name}){
      	croak "ERROR: Trying to merge contigs with conflicting data:\n".Dumper($contigs->{$name})."\n".Dumper($contig)
          unless $contig->compare($contigs->{$name});
      } else {
      	$contigs->{$name} = $contig;
      }
    }
  }
  return $contigs;
}

sub parse_samples{
  my($header_txt,$study,$protocol,$accession,$accession_source,$description) = @_;

  my $samples = {};

  foreach my $line (split(/\n/,$header_txt)){
    if($line =~ /^\@RG/){

      my ($name) = $line =~ /SM:([^\t]+)/;
      my ($platform) = $line =~ /PL:([^\t]+)/;

      $samples->{$name} = new Sanger::CGP::Vcf::Sample(
        -name => $name,
        -study => $study,
        -platform => $platform,
        -seq_protocol => $protocol,
        -accession => $accession,
        -accession_source => $accession_source,
        -description => $description
      ) unless exists $samples->{$name};
    }
  }
  return $samples;
}

1;
