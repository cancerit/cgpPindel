#!/usr/bin/perl

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

my $config_path;
BEGIN {
	use Cwd qw(abs_path);
	use File::Basename;
	my $prog_path = abs_path($0);
	unshift @INC, dirname($prog_path).'/../lib';
};
use strict;
use English qw( -no_match_vars );
use warnings;
use Carp;
use Getopt::Long;
use Pod::Usage;
use Cwd qw(abs_path);
use File::Basename qw(fileparse);
use File::Temp qw(tempfile);
use File::Spec::Functions qw(splitpath);
use Vcf;
use Tabix;

use Sanger::CGP::PindelPostProcessing::VcfSoftFlagger;
use Sanger::CGP::PindelPostProcessing::AbstractExe qw(get_version);

eval {
  main(option_builder(@ARGV));
  1;
} or do {
  warn $EVAL_ERROR if($EVAL_ERROR);
  warn $CHILD_ERROR if($CHILD_ERROR);
  warn $OS_ERROR if($OS_ERROR);
  croak "Something went wrong, check messages\n";
};

sub option_builder{
	my @args = @_;
	my %opts = ();
	my @random_args = ();
	my $result = &GetOptions (
		'h|help' => \$opts{'h'},
		'i|input=s' => \$opts{'i'},
		'o|output=s' => \$opts{'o'},
		'r|rules=s' => \$opts{'r'},
		'sr|softrules=s' => \$opts{'sr'},
		'a|annot=s' => \$opts{'a'},
		'u|unmatched=s' => \$opts{'u'},
		's|simrep=s' => \$opts{'s'},
		'p|processid=s' => \$opts{'p'},
		'v|version' => \$opts{'v'},
		'<>' => sub{push(@random_args,shift(@_));}
	);

	pod2usage(0) if $opts{'h'};

  if($opts{'v'}) {
    print 'Version: ',Sanger::CGP::PindelPostProcessing::VcfSoftFlagger->VERSION,"\n";
    exit 0;
  }

	pod2usage("Unrecognised command line arguments: ".join(",",@random_args)) if(scalar(@random_args));

	pod2usage(q{'-r' has not been defined}) unless(defined $opts{'r'});
	pod2usage(q{'-r' is not a valid file}) unless(-f $opts{'r'});
	pod2usage(q{'-a' has not been defined}) unless(defined $opts{'a'});
	pod2usage(q{'-a' is not a valid file}) unless(-f $opts{'a'});
	pod2usage(q{'-u' has not been defined}) unless(defined $opts{'u'});
	pod2usage(q{'-u' is not a valid file}) unless(-f $opts{'u'});
	pod2usage(q{'-s' has not been defined}) unless(defined $opts{'s'});
	pod2usage(q{'-s' is not a valid file}) unless(-f $opts{'s'});

	##optional
	pod2usage(q{'-sr' is not a valid file}) if($opts{'sr'} && ! -f $opts{'sr'});

	return \%opts;
}

sub main{
	my ($opts) = @_;
	#validate_file($opts->{i}); #validate the incomming file

	$ENV{VCF_FLAGGING_UNMATCHED_NORMALS} = $opts->{'u'};
	$ENV{VCF_FLAGGING_REPEATS} = $opts->{'s'};
	$ENV{VCF_IS_CODING} = $opts->{'a'};

	my $flagger = new Sanger::CGP::PindelPostProcessing::VcfSoftFlagger(
		filter_rules => $opts->{'r'},
		flag_rules => $opts->{'sr'});

	my $vcf = Vcf->new(file => $opts->{i});
	$vcf->parse_header();

	my $out_fh;
	open $out_fh, '>' ,$opts->{'o'} or croak 'Cannot open '.$opts->{'o'}." for writing: $!";

	add_process_log($opts,$vcf);

	$flagger->flag_vcf($vcf,$out_fh);

	close $out_fh;
	$vcf->close();

	validate_file($opts->{o});
}

sub add_process_log {
	my ($opts,$vcf) = @_;

	my $process_version = get_version() || '.';

	$vcf->add_header_line({'key'=>'vcfProcessLog',
			InputVCF => _trim_file_path($opts->{'i'}),
			InputVCFSource => 'FlagVcf.pl',
			InputVCFVer => $process_version,
			InputVCFParam => {
				rules => _trim_file_path($opts->{'r'}),
				annot => _trim_file_path($opts->{'a'}),
				unmatched => _trim_file_path($opts->{'u'}),
				simrep => _trim_file_path($opts->{'s'})},
		},
		'append'=>1);

	$vcf->add_header_line({'key'=>'source', 'value' => 'FlagVcf.pl'}, 'append' => 1);
	$vcf->add_header_line({'key'=>'cgpAnalysisProc', 'value' => $opts->{'p'} }, 'append' => 1) if $opts->{'p'};
}

sub _trim_file_path{
	my ( $string ) = @_;
	my @bits = (split("/", $string));
	return pop @bits;
}

sub validate_file{
	my ($input_file_path) = @_;

	##############################
	## validate the vcf file
	##############################
	eval{
		Vcf::validate($input_file_path);
	};if($@){
		croak($@);
	}
	##############################
}

__END__

=head1 NAME

FlagVcf.pl - Takes a Pindel VCF file and applies a set of flags through a system call to vcf-annotate.

=head1 SYNOPSIS

  General Options:

    --help        (-h)       Brief documentation
    --version     (-v)       Version

  (REQUIRED)
    --rules       (-r)       Full path to a rules file.
    --annot       (-a)       Full path to an indexed tabix annotation file.
    --unmatched   (-u)       Full path to an indexed tabix unmatched normal file.
    --simrep      (-s)       Full path to an indexed tabix simple repeat file.

    --input       (-i)       Full path to input file. Will handle decompression of *.gz automatically using bgzip.
    --output      (-o)       Full path to output file. The output will NOT be bgzipped even if the .gz suffix is set.

  (OPTIONAL)
    --softrules   (-sr)      Full path to a rules file. These are soft filters and as such are assigned to the INFO field.
    --processid   (-p)       The value to appear in cgpAnalysisProc header item. This will overwrite any existing cgpAnalysisProc value.
  Example:

    perl FlagVcf.pl -i mypinfile.vcf -o mypinfile_flagged.vcf -r/VcfFilterRules/PindelFilterRulesPulldown.pm -a /nfs/users/nfs_j/jwh/temp/full_annotation.sorted.gff3.gz -u /nfs/users/nfs_j/jwh/temp/unmatched_normal_pulldown.gff.gz

  Dependances:

    Perl5Lib: Vcf Tabix

    Path Items: bgzip vcf-annotate

=cut
