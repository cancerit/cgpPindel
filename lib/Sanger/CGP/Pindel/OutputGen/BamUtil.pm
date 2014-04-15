package Sanger::CGP::Pindel::OutputGen::BamUtil;

use Sanger::CGP::Pindel::OutputGen;
use Sanger::CGP::Vcf::Contig;
use Sanger::CGP::Vcf::Sample;
our $VERSION = Sanger::CGP::Pindel::OutputGen->VERSION;

use strict;
use Carp;

sub parse_contigs{
  my($header_txt) = @_;
  
  my $contigs = {};
  
  foreach my $line (split(/\n/,$header_txt)){
    if($line =~ /^\@SQ/){
      my ($name) = $line =~ /SN:([^\t]+)/;
      my ($length) = $line =~ /LN:([^\t]+)/;
      my ($assembly) = $line =~ /AS:([^\t]+)/;
      my ($species) = $line =~ /SP:([^\t]+)/;
      
      my $contig = new Sanger::CGP::Vcf::Contig(
        -name => $name,
        -length => $length,
        -assembly => $assembly,
        -species => $species
      );
      
      if(exists $contigs->{$name}){
      	croak "Trying to merge contigs with conflicting data:\n".Dumper($contigs->{$name})."\n".Dumper($contig)
          unless $contigs->{$name}->compare($contig);
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