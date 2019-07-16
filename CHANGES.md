# CHANGES

## 3.2.2

* Handle Input files that may have no reads at all, specifically an issue when generating a normal panel.

## 3.2.1

* Added Dockerfile and docker documentation

## 3.2.0

* Tabix search for high depth/excluded regions now performed in memory using IntervalTrees
  * Reduces runtime of input step by ~50%
  * Improved disk access profile
  * Zero impact on results

## 3.1.2

* 3.0.5 introduced species parsing bug causing single word species names to be invalid.

## 3.1.1

* Fix regression - ability to cope with chromosomes with no events.

## 3.1.0

* Incorporates updated pindel which improves sensitivity
* Internally interpret QCFAIL to determine if whole pair fails

## 3.0.6

* Fixed version tag

## 3.0.5

* Handles species names with spaces in it
* modified checks for species,assembly and checksum

## 3.0.4

* Output bug for pindel BAM/CRAM corrected.  When more than 1 chr in output files had no reads.

## 3.0.3

* Changes to how germline filter determined resulted in dummy germline bed file not being generated as previously.
* This release reinstates the old behaviour.

## 3.0.2

* Correct example rule files for *Fragment.lst files to use FFnnn filter types

## 3.0.1

* Update tabix calls to directly use query_full (solves GRCh38 contig name issues).

## 3.0.0

* Germline bed file is now merged for adjacent regions (#31)
* More compressed intermediate files (#55)
* Change to `Const::Fast` where appropriate (#41)
* Removed TG VG from genotype.
  * Readgroups are always variable, often 1 in data from last few years
  * Not used by our filters.
* Supports BAM/CRAM inputs
* Output will be aligned with inputs
  * bam vs cram
  * bai vs csi
* Although ground work for csi input/output has been done `Bio::DB::HTS` doesn't support csi indexed input yet.
  * Created our own fork at [`cancerit/Bio::DB::HTS`][cancerit-biodbhts] so that this could be enabled.
  * You will need to install this manually or use one of our images for this functionallity.
    * [dockstore-cgpwxs][ds-cgpwxs-git]
    * [dockstore-cgpwxs][ds-cgpwgs-git]

<!-- -->
[cancerit-biodbhts]: https://github.com/cancerit/Bio-DB-HTS/releases/tag/v2.10-rc1
[ds-cgpwxs-git]: https://github.com/cancerit/dockstore-cgpwxs
[ds-cgpwgs-git]: https://github.com/cancerit/dockstore-cgpwgs

## 2.2.5

* Update tabix->query to tabix->query_full

## 2.2.4

* Force sorting of FILTER field to make records easier to diff.
* Fix sorting of final VCF to handle events with same start better when using comparison tools

## 2.2.3

Correct read sorting during collection of DI events.  Caused some events to be split into many and
others to be missed (Thanks to @liangkaiye for patch)

## 2.2.3

Correct read sorting during collection of DI events.  Caused some events to be split into many and
others to be missed (Thanks to @liangkaiye for patch)

## 2.2.2

Correction to sorting of VCF files

## 2.2.0

Reduces the amount of temporary space required and overall I/O

To process 40 million readpairs (40x Tumour + 40x Normal, chr21, 100bp reads):

Original time:

```
User time (seconds): 3553.88
System time (seconds): 63.92
Percent of CPU this job got: 159%
Elapsed (wall clock) time (h:mm:ss or m:ss): 37:51.63
File system inputs: 64
File system outputs: 1782080
```

New time:

```
User time (seconds): 3572.21
System time (seconds): 74.06
Percent of CPU this job got: 167%
Elapsed (wall clock) time (h:mm:ss or m:ss): 36:15.01
File system inputs: 0
File system outputs: 1139128
```

```
Original peak size: 650MB
     New peak size: 291MB
```

__~55%__ reduction in working space and about __40%__ fewer writes to the file system.

Exactly the same results:

```bash
$ diff old/f9c3bc8e-dbc4-1ed0-e040-11ac0d4803a9_vs_f9c3bc8e-dbc1-1ed0-e040-11ac0d4803a9.germline.bed new/f9c3bc8e-dbc4-1ed0-e040-11ac0d4803a9_vs_f9c3bc8e-dbc1-1ed0-e040-11ac0d4803a9.germline.bed

$ diff_bams -a old/f9c3bc8e-dbc4-1ed0-e040-11ac0d4803a9_vs_f9c3bc8e-dbc1-1ed0-e040-11ac0d4803a9_wt.bam -b new/f9c3bc8e-dbc4-1ed0-e040-11ac0d4803a9_vs_f9c3bc8e-dbc1-1ed0-e040-11ac0d4803a9_wt.bam
Reference sequence count passed
Reference sequence order passed
Matching records: 194543

$ diff_bams -a old/f9c3bc8e-dbc4-1ed0-e040-11ac0d4803a9_vs_f9c3bc8e-dbc1-1ed0-e040-11ac0d4803a9_mt.bam -b new/f9c3bc8e-dbc4-1ed0-e040-11ac0d4803a9_vs_f9c3bc8e-dbc1-1ed0-e040-11ac0d4803a9_mt.bam
Reference sequence count passed
Reference sequence order passed
Matching records: 239737

$ /software/CGP/canpipe/live/bin/canpipe_live vcftools --gzvcf old/f9c3bc8e-dbc4-1ed0-e040-11ac0d4803a9_vs_f9c3bc8e-dbc1-1ed0-e040-11ac0d4803a9.flagged.vcf.gz --gzdiff new/f9c3bc8e-dbc4-1ed0-e040-11ac0d4803a9_vs_f9c3bc8e-dbc1-1ed0-e040-11ac0d4803a9.flagged.vcf.gz
...
Comparing individuals in VCF files...
N_combined_individuals:	2
N_individuals_common_to_both_files:	2
N_individuals_unique_to_file1:	0
N_individuals_unique_to_file2:	0
Comparing sites in VCF files...
Found 15321 SNPs common to both files.
Found 0 SNPs only in main file.
Found 0 SNPs only in second file.
After
