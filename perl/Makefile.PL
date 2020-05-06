#!/usr/bin/perl

########## LICENCE ##########
# Copyright (c) 2014-2019 Genome Research Ltd.
#
# Author: CASM/Cancer IT <cgphelp@sanger.ac.uk>
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


use ExtUtils::MakeMaker;

WriteMakefile(
  NAME          => 'cgpPindel',
  LICENSE       => 'agpl_3', # http://search.cpan.org/~dagolden/CPAN-Meta-2.142690/lib/CPAN/Meta/Spec.pm#license
  VERSION_FROM  => 'lib/Sanger/CGP/Pindel.pm',
  EXE_FILES     => [qw( bin/pindel.pl
                        bin/pindel_input_gen.pl
                        bin/pindel_2_combined_vcf.pl
                        bin/FlagVcf.pl
                        bin/pindel_merge_vcf_bam.pl
                        bin/pindel_np_from_vcf.pl
                        bin/pindel_germ_bed.pl
                        bin/pindelCohort.pl
                        bin/pindelCohort_to_vcf.pl
                        bin/pindel_vcfSortNsplit.pl
                        bin/pindel_blat_vaf.pl
                        )],
  PREREQ_PM     => {
                     'Const::Fast' => 0.014,
                     'Try::Tiny' => 0.19,
                     'File::Which' => 0.05,
                     'List::Util' => 1.38,
                     'Bio::DB::HTS' => 1.11,
                     'List::MoreUtils' => 0.33,
                     'Test::Fatal' => 0.013,
                     'Devel::Cover' => 1.09,
                     'Pod::Coverage' => 0.23,
                     'PerlIO::gzip' => 0.20,
                     'Set::IntervalTree' => 0.12,
                   }
);
