LICENCE
=======

Copyright (c) 2014 Genome Research Ltd.

Author: Keiran Raine <cgpit@sanger.ac.uk>

This file is part of cgpPindel.

cgpPindel is free software: you can redistribute it and/or modify it under
the terms of the GNU Affero General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option) any
later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
details.

You should have received a copy of the GNU Affero General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.

cgpPindel
=========

cgpPindel contains the Cancer Genome Projects workflow for [pindel](http://gmt.genome.wustl.edu/pindel/current/).

The is a lightly modified version of pindel v2.0 with CGP specific processing for:

* Input file generation
* Conversion from pindel text output to:
** tumour and normal BAM alignment files
** VCF
** Application of VCF filters.

---

###Dependencies/Install

Please install the following first:

* [PCAP-core](http://github.com/ICGC-TCGA-PanCancer/PCAP-core/releases)
* [cgpVcf](http://github.com/cancerit/cgpVcf/releases)

Please see these for any child dependencies.

Once complete please run:

./setup.sh /some/install/location

Please be aware that this expects basic C compilation libraries and tools to be available,
most are listed in `INSTALL`.

---

##Creating a release
####Preparation
* Commit/push all relevant changes.
* Pull a clean version of the repo and use this for the following steps.

####Cutting the release
1. Update `perl/lib/Sanger/CGP/Pindel.pm` to the correct version (adding rc/beta to end if applicable).
2. Update `Changes` to show major items.
3. Run `./prerelease.sh`
4. Check all tests and coverage reports are acceptable.
5. Commit the updated docs tree and updated module/version.
6. Push commits.
7. Use the GitHub tools to draft a release.
