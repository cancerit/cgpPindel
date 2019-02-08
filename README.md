# cgpPindel

cgpPindel contains the Cancer Genome Projects workflow for [Pindel][pindel-core].

[![Quay Badge][quay-status]][quay-repo]

| Master                                        | Develop                                         |
| --------------------------------------------- | ----------------------------------------------- |
| [![Master Badge][travis-master]][travis-base] | [![Develop Badge][travis-develop]][travis-base] |

The is a lightly modified version of pindel v2.0 with CGP specific processing for:

* Input file generation
* Conversion from pindel text output to:
  * tumour and normal BAM alignment files
  * VCF
  * Application of VCF filters.

Contents:

<!-- TOC depthFrom:2 depthTo:6 withLinks:1 updateOnSave:1 orderedList:0 -->

- [cgpPindel](#cgppindel)
  - [Docker, Singularity and Dockstore](#docker-singularity-and-dockstore)
  - [Dependencies/Install](#dependenciesinstall)
  - [Creating a release](#creating-a-release)
    - [Preparation](#preparation)
    - [Release process](#release-process)
      - [Code changes](#code-changes)
      - [Docker image](#docker-image)
      - [Cutting the release](#cutting-the-release)
  - [LICENCE](#licence)

<!-- /TOC -->

## Docker, Singularity and Dockstore

There are pre-built images containing this codebase on quay.io.

* [cgpPindel][cgpPindel-git]: Contained within this repository - contains the cgpPindel package
* [dockstore-cgpwxs][ds-cgpwxs-git]: Contains tools specific to WXS analysis.
* [dockstore-cgpwgs][ds-cgpwgs-git]: Contains additional tools for WGS analysis.

These were primarily designed for use with dockstore.org but can be used as normal containers.

The docker images know to work correctly after import into a singularity image.

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

### Release process

This project is maintained using HubFlow.

#### Code changes

1. Make appropriate changes
2. Update `perl/lib/Sanger/CGP/Pindel.pm` to the correct version (adding rc/beta to end if applicable).
3. Update `CHANGES.md` to show major items.
4. Run `./prerelease.sh`
5. Check all tests and coverage reports are acceptable.
6. Commit the updated docs and updated module/version.
7. Push commits.

#### Docker image

1. Use the GitHub tools to draft a release.
2. Build image locally
3. Run example inputs and verify any changes are acceptable
4. Bump version in `Dockerfile`
5. Push changes

#### Cutting the release

1. Check state on Travis
2. Generate the release (add notes to GitHub)
3. Confirm that image has been built on [quay.io][quay-builds]
4. Update the [dockstore][dockstore-cgpPindel] entry, see [their docs][dockstore-get-started].

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
[cgpPindel-git]: https://github.com/cancerit/cgpPindel
[pindel-core]: http://gmt.genome.wustl.edu/pindel/current

<!-- Travis -->
[travis-base]: https://travis-ci.org/cancerit/cgpPindel
[travis-master]: https://travis-ci.org/cancerit/cgpPindel.svg?branch=master
[travis-develop]: https://travis-ci.org/cancerit/cgpPindel.svg?branch=dev

<!-- Quay.io -->
[quay-status]: https://quay.io/repository/wtsicgp/cgpPindel/status
[quay-repo]: https://quay.io/repository/wtsicgp/cgpPindel
[quay-builds]: https://quay.io/repository/wtsicgp/cgpPindel?tab=builds

<!-- Dockstore -->
[dockstore-cgpPindel]: https://dockstore.org/containers/quay.io/wtsicgp/cgpPindel
