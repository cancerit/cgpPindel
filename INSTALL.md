# Installation

    ./setup.sh /path/to/installation

`/path/to/installation` is where you want the `bin`, `lib` folders to be created. The following tools will be installed: `cgpPindel` and `pindel`.

⚠️ *This distribution will only works on `*NIX` type systems.*

# System Dependencies

**Perl:** Minimum version is `5.10.1` (tested with `5.16.3`).

<!-- we should not duplicate this info -->
* Ubuntu 16.04: System dependencies should be satisfied with a successful installation of `PCAP-core` and `cgpVcf`. If you find otherwise please let us know! See [`Dockerfile`](Dockerfile).

* CentOS 6.4:

        yum install \
            zlib-devel \
            gcc-c++ \
            ncurses-devel.x86_64
