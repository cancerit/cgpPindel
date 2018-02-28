# cgpPindel

cgpPindel contains the Cancer Genome Projects workflow for [Pindel][pindel-core].

| Master                                        | Develop                                         |
| --------------------------------------------- | ----------------------------------------------- |
| [![Master Badge][travis-master]][travis-base] | [![Develop Badge][travis-develop]][travis-base] |

The is a lightly modified version of pindel v2.0 with CGP specific processing for:

* Input file generation
* Conversion from pindel text output to:
    * tumour and normal BAM alignment files
    * VCF
    * Application of VCF filters.

## Docker, Singularity and Dockstore

There are pre-built images containing this codebase on quay.io.

* [dockstore-cgpwxs][ds-cgpwxs-git]
  * Contains tools specific to WXS analysis.
* [dockstore-cgpwgs][ds-cgpwgs-git]
  * Contains additional tools for WGS analysis.

These were primarily designed for use with dockstore.org but can be used as normal containers.

The docker images are know to work correctly after import into a singularity image.

## Dependencies/Install

Please install the following first:

* [PCAP-core v2.0+][pcap-core-rel]
* [cgpVcf v2.0+][cgpvcf-rel]

Please see these for any child dependencies.

Once complete please run:

```
./setup.sh /some/install/location
```

Please use `setup.sh` to install any other dependencies.  Setting the environment variable
`CGP_PERLLIBS` allows you to to append to `PERL5LIB` during install.  Without this all dependancies
are installed into the target area.

Please be aware that this expects basic C compilation libraries and tools to be available.

## Creating a release

### Preparation

* Commit/push all relevant changes.
* Pull a clean version of the repo and use this for the following steps.

###Cutting the release

1. Update `perl/lib/Sanger/CGP/Pindel.pm` to the correct version (adding rc/beta to end if applicable).
1. Update `CHANGES.md` to show major items.
1. Run `./prerelease.sh`
1. Check all tests and coverage reports are acceptable.
1. Commit the updated docs and updated module/version.
1. Push commits.
1. Use the GitHub tools to draft a release.

## LICENCE

```
Copyright (c) 2014-2018 Genome Research Ltd.

Author: CASM/Cancer IT <cgphelp@sanger.ac.uk>

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

1. The usage of a range of years within a copyright statement contained within
this distribution should be interpreted as being equivalent to a list of years
including the first and last year specified and all consecutive years between
them. For example, a copyright statement that reads ‘Copyright (c) 2005, 2007-
2009, 2011-2012’ should be interpreted as being identical to a statement that
reads ‘Copyright (c) 2005, 2007, 2008, 2009, 2011, 2012’ and a copyright
statement that reads ‘Copyright (c) 2005-2012’ should be interpreted as being
identical to a statement that reads ‘Copyright (c) 2005, 2006, 2007, 2008,
2009, 2010, 2011, 2012’."
```

<!-- References -->
[cgpvcf-rel]: https://github.com/cancerit/cgpVcf/releases
[pcap-core-rel]: https://github.com/cancerit/PCAP-core/releases
[ds-cgpwxs-git]: https://github.com/cancerit/dockstore-cgpwxs
[ds-cgpwgs-git]: https://github.com/cancerit/dockstore-cgpwgs
[pindel-core]: http://gmt.genome.wustl.edu/pindel/current

<!-- Travis -->
[travis-base]: https://travis-ci.org/cancerit/cgpPindel
[travis-master]: https://travis-ci.org/cancerit/cgpPindel.svg?branch=master
[travis-develop]: https://travis-ci.org/cancerit/cgpPindel.svg?branch=dev
