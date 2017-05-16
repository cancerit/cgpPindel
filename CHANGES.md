### 2.2.2
Correction to sorting of VCF files

### 2.2.0
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

__~55% reduction in working space and about 40% fewer writes to the file system.__

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
After filtering, kept 16309 out of a possible 16309 Sites
Run Time = 6.00 seconds
```

### 2.0.4
* Permits empty results files

### 2.0.0
* Migrates all Tabix and Bio::DB::Sam to Bio::DB::HTS::Tabix and Bio::DB::Sam
* Cleans up install
