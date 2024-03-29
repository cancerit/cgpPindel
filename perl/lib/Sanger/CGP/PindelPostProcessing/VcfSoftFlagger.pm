# Copyright (c) 2014-2021 Genome Research Ltd
#
# Author: CASM/Cancer IT <cgphelp@sanger.ac.uk>
#
# This file is part of cgpPindel.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
# 1. The usage of a range of years within a copyright statement contained within
# this distribution should be interpreted as being equivalent to a list of years
# including the first and last year specified and all consecutive years between
# them. For example, a copyright statement that reads ‘Copyright (c) 2005, 2007-
# 2009, 2011-2012’ should be interpreted as being identical to a statement that
# reads ‘Copyright (c) 2005, 2007, 2008, 2009, 2011, 2012’ and a copyright
# statement that reads ‘Copyright (c) 2005-2012’ should be interpreted as being
# identical to a statement that reads ‘Copyright (c) 2005, 2006, 2007, 2008,
# 2009, 2010, 2011, 2012’.
#
package Sanger::CGP::PindelPostProcessing::VcfSoftFlagger;

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

	## We also need to update the vcf process log to represent the addition of flagging...
	## currently this wil have to be done outside of VcfAnnotate as vcfProcessLog is not a standard header (TCGA).

	my $filter_file = $args{filter_rules};
	my $flag_file = $args{flag_rules};
	my $zcat_exe = $args{zcat_exe} || 'zcat';
	my $bgzip_exe = $args{bgzip_exe} || 'bgzip';

    my $self = {
    	filter_file => $filter_file,
    	flag_file => $flag_file,
    	zcat_exe => $zcat_exe,
    	bgzip_exe => $bgzip_exe
    };
    bless $self, $class;
    return $self;
}

=head3 flag

=over

=item Param 1 full input vcf file path. If the path ends in '.gz' an attempt to un-bgzip the file before processing will be made.

=item Param 2 full output vcf file path. If the path ends in '.gz' an attempt to bgzip the file after processing will be made.

=back

=cut
sub flag{
	my($self,$input,$output) = @_;

	## build the comand line...
	croak("Input |$input| is not a valid file") unless (-f $input);
	croak("Output has not been defined") unless ( $output );

	croak("No filter rule file defined") unless ( $self->{filter_file} );
	croak("Not a valid file path: ".$self->{filter_file}) unless ( -f $self->{filter_file} );

	open(my $out_fh, ">", $output) || croak("could not open $output for writing: $!");

	my $vcf = Vcf->new(file => $input);
	$vcf->parse_header();

	$self->flag_vcf($vcf,$out_fh);

	close $out_fh;
	$vcf->close();
}

=head3 flag

=over

=item Param 1 Vcf object.

=item Param 2 an open writable file handle.

=back

=cut
sub flag_vcf{
	my ($self,$vcf,$out_fh) = @_;

	## read in the rules to the soft flagger. Doe this for both filters and flags.
	$self->_parse_user_defined_filters($vcf,$self->{'filter_file'},'udef_filters') if(defined $self->{'filter_file'});
	$self->_parse_user_defined_filters($vcf,$self->{'flag_file'},'udef_flags') if(defined $self->{'flag_file'});

	## print out the header.
	print $out_fh reformat_header($vcf);

 	## print out the body of the file.
	while (my $line = $vcf->next_line()){
		chomp $line;
		my $line_bits = [split(/\t/,$line)];
		my $new_rec = $self->apply_user_defined_filters($vcf,$line_bits,'udef_filters');
		$new_rec = $self->apply_user_defined_filters($vcf,$new_rec,'udef_flags');
		print $out_fh $vcf->format_line($new_rec);
	}
}

=head3 reformat_header

  Vcf module makes a mess of this and does not group common header types
  simply take the existing data and use the order that each header type is first seen
  and build an array of each following element of that type

=over

=item Param 1 Vcf object.

=item Returns the header as a formatted string.

=back

=cut
sub reformat_header {
	my $vcf = shift;
	my @lines = split /\n/, $vcf->format_header;
	my $format_head = shift @lines;
	my $col_head = pop @lines;
	my %header_sets;
	for my $h_line(@lines) {
		my ($h_type) = $h_line =~ m/^\#\#([^=]+)/;
		unless(exists $header_sets{$h_type}) {
			$header_sets{$h_type} = [];
		}
		push @{$header_sets{$h_type}}, $h_line;
	}

	my @final_header;
	push @final_header, $format_head;
	for my $h_type(sort keys %header_sets) {
		push @final_header, @{$header_sets{$h_type}};
	}
	push @final_header, $col_head;
	return (join "\n", @final_header)."\n";
}

=head3 _parse_user_defined_filters

  Modified from vcf-annotate vcftools. Takes a vcf-annotate filter file and reads it into the flagger.

=over

=item Param 1 Vcf object.

=item Param 2 Full file path to a filter file.

=item Param 3 Filter type string (udef_filters | udef_flags). This determins how the filter file is treated i.e. as INFO flags or filters.

=back

=cut
sub _parse_user_defined_filters
{
	my ($self,$vcf,$file,$root_filter_set_name) = @_;

  my $filters = [];
  open my $IN, '<', $file or croak "Failed to read $file: $!";
  my $module = <$IN>;
  chomp $module;
  eval "require $module";
  while(my $flag = <$IN>) {
    chomp $flag;
    push @{$filters}, $module->rule($flag);
  }
  close $IN;

	my $info_descs = [];

	delete($$self{$root_filter_set_name}); ## reset the filterlist.

	for my $filter (@$filters){

		if ( !exists($$filter{tag}) ) { croak("Missing 'tag' key for one of the filters in $file\n"); }
		if ( $$filter{tag}=~m{^INFO/(.+)$} ) { $$filter{info_tag} = $1; }
		elsif ( $$filter{tag}=~m{^FORMAT/(.+)$} ) { $$filter{format_tag} = $1; }
		else { croak("Currently only INFO and FORMAT tags are supported. Could not parse the tag [$$filter{tag}]\n"); }

		if ( !exists($$filter{name}) ) { croak("Missing 'name' key for the filter [$$filter{tag}]\n"); }
		if ( !exists($$filter{desc}) ) { croak("Missing 'desc' key for the filter [$$filter{tag}]\n"); }

		my $name = $$filter{name};

		if($root_filter_set_name eq 'udef_flags'){
			$vcf->add_header_line({key=>'INFO',ID=>$name,Number=>0,Type=>'Flag',Description=>$$filter{desc}},silent=>1);
		}else{
			$vcf->add_header_line({key=>'FILTER',ID=>$name,Description=>$$filter{desc}},silent=>1);
		}

		if ( !exists($$filter{apply_to}) or lc($$filter{apply_to}) eq 'all' ){
			$$self{$root_filter_set_name}{'all'}{$name} = $filter;
			$$self{$root_filter_set_name}{'s'}{$name}   = $filter;
			$$self{$root_filter_set_name}{'i'}{$name}   = $filter;
		}elsif ( exists($$filter{apply_to}) and lc($$filter{apply_to}) eq 'snps' ){
			$$self{$root_filter_set_name}{'s'}{$name}   = $filter;
			$$self{udef_filters_typecheck_needed} = 1;
		}elsif ( exists($$filter{apply_to}) and lc($$filter{apply_to}) eq 'indels' ){
			$$self{$root_filter_set_name}{'i'}{$name}   = $filter;
			$$self{udef_filters_typecheck_needed} = 1;
		}
	}
}

=head3 apply_user_defined_filters

  Modified from vcf-annotate vcftools. Applies the filters of a specified filter-set to a line of data.

=over

=item Param 1 Vcf object.

=item Param 2 Array-ref representing the vcf line being worked on.

=item Param 3 Filter type string (udef_filters | udef_flags). This determins which filter-set is to be applied i.e. as INFO flags or filters.

=item Returns An array-ref of the modified line.

=back

=cut
sub apply_user_defined_filters{
	my ($self,$vcf,$line,$root_filter_set_name) = @_;

	my($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF);
	$CHROM  = $$line[0];
	$POS	= $$line[1];
	$FAIL   = 0;
	$PASS   = 1;
	$RECORD = $line;
	$VCF	= $vcf;

	my %filters = ();
	my %apply = ();

	if(exists($$self{$root_filter_set_name})){
		if ( $$self{udef_filters_typecheck_needed} ){
			# Check if the line has an indel, SNP or both
			for my $alt (split(/,/,$$line[4])){
				my ($type,$len,$ht) = $vcf->event_type($$line[3],$alt);
				if ( exists($$self{$root_filter_set_name}{$type}) )
				{
					%filters = ( %filters, %{$$self{$root_filter_set_name}{$type}} );
				}
			}
			# Return if the line does not have the wanted variant type
			if ( !scalar %filters ) { return $line; }
		}else{
			%filters = %{$$self{$root_filter_set_name}{all}};
		}

		for my $filter (values %filters){
			if ( exists($$filter{info_tag}) ){
				$MATCH = $vcf->get_info_field($$line[7],$$filter{info_tag});
				if ( !defined $MATCH ) { next; }
			}elsif ( exists($$filter{format_tag}) ){
				my $idx = $vcf->get_tag_index($$line[8],$$filter{format_tag},':');
				if ( $idx < 1 ) { next; }
				$MATCH = $vcf->get_sample_field($line,$idx);
			}

			my $result = &{$$filter{test}}($MATCH,$CHROM,$POS,$FAIL,$PASS,$RECORD,$VCF);

			if($root_filter_set_name eq 'udef_flags'){
				if($result == $PASS){
					$apply{ $$filter{name} } = '';
				}
			}else{
				$apply{ $$filter{name} } = 1 if($result != $PASS);
			}
		}
	}

	if($root_filter_set_name eq 'udef_flags'){
		$$line[7] = $vcf->add_info_field($$line[7],%apply) if(scalar keys %apply);
	}else{
		if(scalar keys %apply) {
			$$line[6] = join ';', sort keys %apply;
		}
		else {
			$$line[6] = 'PASS'
		}
	}

	return $line;
}

sub filter_file{
	my($self,$filter_file) = @_;
	if($filter_file){
		$self->{'filter_file'} = $filter_file;
	}
	return $filter_file;
}

sub flag_file{
	my($self,$flag_file) = @_;
	if($flag_file){
		$self->{'flag_file'} = $flag_file;
	}
	return $flag_file;
}
