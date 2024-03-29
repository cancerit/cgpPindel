# cgpPindel

cgpPindel contains the Cancer Genome Projects workflow for [Pindel][pindel-core].

| Master                                        | Develop                                         |
| --------------------------------------------- | ----------------------------------------------- |
| [![Master Badge][circle-master]][circle-base] | [![Develop Badge][circle-develop]][circle-base] |

[![pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white)](https://github.com/pre-commit/pre-commit)

The is a lightly modified version of pindel v2.0 with CGP specific processing for:

- Input file generation
- Conversion from pindel text output to:
  - tumour and normal BAM alignment files
  - VCF
  - Application of VCF filters.

Details of execution and referencing can be found in the [wiki][cgppindel-wiki]

Contents:

- [Docker, Singularity and Dockstore](#docker-singularity-and-dockstore)
- [Dependencies/Install](#dependenciesinstall)
- [Nextflow](#nextflow)
  - [Setup personal nextflow](#setup-personal-nextflow)
  - [Profiles](#profiles)
  - [Workflow entry points](#workflow-entry-points)
    - [`pindel_pl`](#pindel_pl)
    - [`np_generation`](#np_generation)
  - [Sub-workflows](#sub-workflows)
- [Developers](#developers)
  - [Updating licence headers](#updating-licence-headers)
  - [Code changes](#code-changes)
  - [Testing](#testing)
    - [Regression CI](#regression-ci)
    - [Public CI](#public-ci)
    - [Cutting the release](#cutting-the-release)
- [LICENCE](#licence)

## Docker, Singularity and Dockstore

There are pre-built images containing this codebase on quay.io.  When pulling an image you must specify
the version there is no `latest`.

- [cgpPindel quay.io][quay-repo]: Contained within this repository
  - Smallest build required to use cgpPindel
  - Not linked to Dockstore (yet)
  - Updated most frequently
- [dockstore-cgpwxs][ds-cgpwxs-git]: Contains tools specific to WXS analysis.
- [dockstore-cgpwgs][ds-cgpwgs-git]: Contains additional tools for WGS analysis.

These were primarily designed for use with dockstore.org but can be used as normal containers.

The docker images are known to work correctly after import into a singularity image.

## Dependencies/Install

When doing a native install please install the following first:

- [PCAP-core v2.0+][pcap-core-rel]
- [cgpVcf v2.0+][cgpvcf-rel]

Please see these for any child dependencies.

Once complete please run:

```
./setup.sh /some/install/location
```

Please use `setup.sh` to install any other dependencies.  Setting the environment variable
`CGP_PERLLIBS` allows you to to append to `PERL5LIB` during install.  Without this all dependancies
are installed into the target area.

Please be aware that this expects basic C compilation libraries and tools to be available.

## Nextflow

Initial Nextflow bindings for cgpPindel.

### Setup personal nextflow

If you don't have a central nextflow install this will get you running with a limited environment:

```bash
# seems silly but a python venv is a nice way to handle this during dev
python3 -m venv .venv
source .venv/bin/activate
# compute head nodes may need you to limit Java accessing all memory
export NXF_OPTS="-Xms500M -Xmx2G"
curl get.nextflow.io | bash
mv nextflow .venv/bin/.
```

If you have any issues installing refer to Nextflow documentation, not the issue tracker for this repo.

### Profiles

Refer to nextflow for an explanation of profiles.  The following are available:

- Job management
  - local
    - spawned jobs use execution host
  - lsf
    - spawned jobs are submitted via `bsub`
- Execution method
  - `<none>`
    - Expects to find programs in `PATH`
  - singularity
    - Provide image file via `-with-singularity [singularity image]`
  - docker
    - Provide image via `-with-docker [docker image]`

For example to run on a LSF farm with singularity the profile would be:

```
... -profile lsf,singularity ...
```

While native install with lsf would be:

```
... -profile lsf ...
```

### Workflow entry points

There are 2 top level entry points:

```
... -entry pindel_pl ...
# or
... -entry np_generation ...
```

#### `pindel_pl`

Executes a tumour/normal paired analysis.

CPU and memory are controlled via `nextflow.config`, configs are additive, see nextflow documentation.

For the workflow options run:

```
nextflow -entry pindel_pl --help
```

#### `np_generation`

Executes pindel with a "dummy" tumour for a file listing of input BAMs.

CPU and memory are controlled via `nextflow.config`, configs are additive, see nextflow documentation.

For the workflow options run:

```
nextflow -entry np_generation --help
```

### Sub-workflows

The nextflow code has been implemented with DSL2 so that the workflows can be composed into larger components.

The items above can be addressed for this purpose via:

- `subwf_pindel_pl`
- `subwf_np_gen`

## Developers

Please use `pre-commit` on this project.  You can install to `$HOME/bin` via:

```bash
curl https://pre-commit.com/install-local.py | python -
```

In you checkout please run:

```bash
pre-commit install
```

### Updating licence headers

Please use [skywalking-eyes](https://github.com/apache/skywalking-eyes).

Expected workflow:

```bash
# recent build, change to apache/skywalking-eyes:0.2.0 once released
export DOCKER_IMG=ghcr.io/apache/skywalking-eyes/license-eye
```

1. Check state before modifying `.licenserc.yaml`:
   - `docker run -it --rm -v $(pwd):/github/workspace $DOCKER_IMG header check`
   - You should get some 'valid' here, those without a header as 'invalid'
1. Modify `.licenserc.yaml`
1. Apply the changes:
   - `docker run -it --rm -v $(pwd):/github/workspace $DOCKER_IMG header fix`
1. Add/commit changes

This is executed in the CI pipeline.

*DO NOT* edit the header in the files, please modify the date component of `content` in `.licenserc.yaml`.  The only exception being:

- `README.md`

If you need to make more extensive changes to the license carefully test the pattern is functional.

### Code changes

This project is maintained using the [HubFlow][hubflow-docs] methodology.

1. Make appropriate changes
1. Update `perl/lib/Sanger/CGP/Pindel.pm` to the correct version (adding rc/beta to end if applicable).
1. Update `CHANGES.md` to show major items.
1. Commit the updated docs and updated module/version.
1. Push commits.

### Testing

#### Regression CI

An internal CI system is used to validate each release using real, large scale datasets.

#### Public CI

Circleci is used to:

- Build Docker image (unit tests are part of build)
- Validate expected tools exist
- For tags only: push image to quay.io

CI only runs for:

- Branches with pull-requests
- Default branch (`dev`)
- Tags

#### Cutting the release

Internal regression CI processes must be completed prior to this.

1. Check state on [Circleci][circle-repo]
1. Generate the release (add notes to GitHub)
1. Confirm that image has been pushed to [quay.io][quay-tags]

## LICENCE

```
Copyright (c) 2014-2021 Genome Research Ltd.

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

<!-- Circle-ci -->

<!-- Quay.io -->

[cgppindel-wiki]: https://github.com/cancerit/cgpPindel/wiki
[cgpvcf-rel]: https://github.com/cancerit/cgpVcf/releases
[circle-base]: https://circleci.com/gh/cancerit/cgpPindel.svg?style=shield
[circle-develop]: https://circleci.com/gh/cancerit/cgpPindel.svg?style=shield&branch=dev%3B
[circle-master]: https://circleci.com/gh/cancerit/cgpPindel.svg?style=shield&branch=master%3B
[circle-repo]: https://app.circleci.com/pipelines/github/cancerit/cgpPindel
[ds-cgpwgs-git]: https://github.com/cancerit/dockstore-cgpwgs
[ds-cgpwxs-git]: https://github.com/cancerit/dockstore-cgpwxs
[hubflow-docs]: https://datasift.github.io/gitflow/
[pcap-core-rel]: https://github.com/cancerit/PCAP-core/releases
[pindel-core]: http://gmt.genome.wustl.edu/pindel/current
[quay-repo]: https://quay.io/repository/wtsicgp/cgppindel
[quay-tags]: https://quay.io/repository/wtsicgp/cgppindel?tab=tags
