# cgpPindel

cgpPindel contains the Cancer Genome Projects workflow for [pindel](http://gmt.genome.wustl.edu/pindel/current/).

| Master                                        | Develop                                         |
| --------------------------------------------- | ----------------------------------------------- |
| [![Master Badge][travis-master]][travis-base] | [![Develop Badge][travis-develop]][travis-base] |

The is a lightly modified version of pindel v2.0 with CGP specific processing for:

* Input file generation
* Conversion from pindel text output to:
    - tumour and normal BAM alignment files
    - VCF
    - Application of VCF filters.

# Contents

- [cgpPindel](#cgppindel)
- [Contents](#contents)
- [Installation](#installation)
- [Contributing](#contributing)
- [License](#license)

# Installation

Please install the following first:

* [PCAP-core v2.0+](http://github.com/cancerit/PCAP-core/releases)
* [cgpVcf v2.0+](http://github.com/cancerit/cgpVcf/releases)

Once complete please run:

    ./setup.sh /some/install/location

⚠️ `cgpPindel` system dependencies should be satisfied with a successful installation of `PCAP-core` and `cgpVcf`. If you find otherwise please let us know!

Setting the environment variable `CGP_PERLLIBS` allows you to to append to `PERL5LIB` during install. Without this all dependancies are installed into the target area.

# Contributing

Contributions are welcome, and they are greatly appreciated, check our [contributing guidelines](CONTROBUTING.md)!

# License

```
Copyright (c) 2014-2016 Genome Research Ltd.

Author: Cancer Genome Project <cgpit@sanger.ac.uk>

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

<!-- Travis -->
[travis-base]: https://travis-ci.org/cancerit/cgpVcf
[travis-master]: https://travis-ci.org/cancerit/cgpVcf.svg?branch=master
[travis-develop]: https://travis-ci.org/cancerit/cgpVcf.svg?branch=develop
