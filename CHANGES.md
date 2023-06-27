# CHANGES

## 1.0.1
- Added a file check after blat step

## 1.0.0
- Added filters to FlagVcf.pl to allow flagging of per-sample vcf outputs
- Fixed bugs in Implement.pm and pindelCohortVafSliceFill.pl
- Adds code to allow single sample processing with more accurate VAF calculations (via BLAT)
- Status of new scripts, "pre-release" indicates defaults and CLI may change:
  - stable
    - pindelCohort.pl
    - pindel_blat_vaf.pl
  - pre-release
    - pindelCohort_to_vcf.pl
    - pindel_vcfSortNsplit.pl
    - pindelCohortMerge.pl
    - pindelCohortVafFill.pl
    - pindelCohortVafSplit.pl
    - pindelCohortVafSliceFill.pl
- pinning to pindel v3.6.0
- Switch license management to skywalking-eyes.
