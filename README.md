cgpVcf
======

cgpVcf contains a set of common perl utilities for generating consistent Vcf headers and an
installation method which resolves some upatched issues in the Tabix package (reported).

It primarily exists to prevent code duplication between some other projects.

---

###Dependencies/Install
Some of the code included in this package has dependencies on several C packages:

 * [Tabix](https://github.com/samtools/tabix)
 * [vcftools](http://vcftools.sourceforge.net/)

And various perl modules.

Please use `setup.sh` to install the dependencies.  Please be aware that this expects basic C
compilation libraries and tools to be available, most are listed in `INSTALL`.

---

##Creating a release
####Preparation
* Commit/push all relevant changes.
* Pull a clean version of the repo and use this for the following steps.

####Cutting the release
1. Update `lib/Sanger/CGP/Vcf.pm` to the correct version (adding rc/beta to end if applicable).
2. Update `Changes` to show major items.
3. Run `./prerelease.sh`
4. Check all tests and coverage reports are acceptable.
5. Commit the updated docs tree and updated module/version.
6. Push commits.
7. Use the GitHub tools to draft a release.
