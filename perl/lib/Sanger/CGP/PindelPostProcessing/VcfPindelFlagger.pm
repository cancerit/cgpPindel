####################################################
# Copyright (c) 2013 Genome Research Ltd.
# Author: Cancer Genome Project, cgpit@sanger.ac.uk
# See LICENCE.TXT for details
####################################################
package Sanger::CGP::PindelPostProcessing::VcfPindelFlagger;

# Author: jwh (Original)

### SVN INFO ###
# $LastChangedDate: $
# $LastChangedRevision: $ # not the tagged version
# $LastChangedBy: $
# $Id: VcfPindelFlagger.pm $
################

use strict;
use Carp;
use English qw( -no_match_vars );
use Sanger::CGP::Pindel;

1;

sub new{
	my $proto = shift;
	my (%args) = @_;
	my $class = ref($proto) || $proto;

	## Need to set the following:
	## 1) The filters file. There will be one of these for each data type. This is not
	## ideal as many of the flags will be repeated across the files. Need to ensure
	## that consistent testing is performed on each file.
	## 2) The unmatched normal panel bed file.
	## 3) The simple repeats bed file.
	## 4) The microsatelite bed file.
	## 5) The annotation bed/gff file.

	## We also need to update the vcf process log to represent the addition of flagging...
	## currently this wil have to be done outside of VcfAnnotate as vcfProcessLog is not a standard header (TCGA).


	my $filter_file = $args{filter_rules};
	my $unmatched_normal_file = $args{unmatched_normals};
	my $simple_repeats_file = $args{simple_repeats};
	my $microsatelites_file = $args{microsatelites};
	my $coding_annotation_file = $args{coding_annotation};
	my $vcf_annotate_exe = $args{vcf_annotate_exe} || 'vcf-annotate';
	my $zcat_exe = $args{zcat_exe} || 'zcat';
	my $bgzip_exe = $args{bgzip_exe} || 'bgzip';

    my $self = {
    	filter_file => $filter_file,
    	unmatched_normal_file => $unmatched_normal_file,
    	simple_repeats_file => $simple_repeats_file,
    	microsatelites_file => $microsatelites_file,
    	coding_annotation_file => $coding_annotation_file,
    	vcf_annotate_exe => $vcf_annotate_exe,
    	zcat_exe => $zcat_exe,
    	bgzip_exe => $bgzip_exe,
    };
    bless $self, $class;
    return $self;
}

=head flag

Given input and output file paths will make a system call to vcf-annotate

=item Param 1 full input vcf file path. If the path ends in '.gz' an attempt to un-bgzip the file before processing will be made.

=item Param 2 full output vcf file path. If the path ends in '.gz' an attempt to bgzip the file after processing will be made.

=cut
sub flag{
	my($self,$input,$output) = @_;

	## build the comand line...
	croak("Input |$input| is not a valid file") unless (-f $input);
	croak("Output has not been defined") unless ( $output );

	croak("No filter rule file defined") unless ( $self->{filter_file} );
	croak("Not a valid file path: ".$self->{filter_file}) unless ( -f $self->{filter_file} );

	my $cmd = "cat  $input | ";
	my $zipped = 0;
	$zipped = 1 if($input =~ /.+\.b?gz$/);
	## cant get this to work..
	#my $out_zipped = 0;
	#$out_zipped = 1 if($output =~ /.+\.gz$/);

	if($zipped){
		$cmd .= $self->{bgzip_exe} .  ' -d | ';
	}

	$cmd .= $self->{vcf_annotate_exe} . " --filter " . $self->{filter_file};
	#$cmd .= " | " . $self->{bgzip_exe} . ' |' if($out_zipped);
	$cmd .= " > $output";

	## set up environment variables to be used by the flagging rules in vcf-annotate.
	$ENV{VCF_FLAGGING_UNMATCHED_NORMALS} = $self->{unmatched_normal_file};
	$ENV{VCF_FLAGGING_REPEATS} = $self->{simple_repeats_file};
	$ENV{VCF_FLAGGING_MICROSATELITE} = $self->{microsatelites_file};
	$ENV{VCF_IS_CODING} = $self->{coding_annotation_file};

	## execute the comand.
	#warn $cmd;
	my $ecode = system($cmd);
	croak("Problem executing command: $cmd") if($ecode);
}
