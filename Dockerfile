FROM quay.io/wtsicgp/pcap-core:5.6.1 as builder

# hadolint ignore=DL3002
USER  root

# ALL tool versions used by opt-build.sh
# need to keep in sync with setup.sh
ENV VER_CGPVCF="v2.2.1"\
    VER_VCFTOOLS="0.1.16"\
    VER_BLAT="v385"

# hadolint ignore=DL3008
RUN apt-get -yq update \
&& apt-get install -yq --no-install-recommends locales g++ make gcc pkg-config zlib1g-dev \
&& locale-gen en_US.UTF-8 \
&& update-locale LANG=en_US.UTF-8

ENV OPT /opt/wtsi-cgp
ENV PATH=$OPT/bin:$OPT/biobambam2/bin:$PATH \
    PERL5LIB=$OPT/lib/perl5 \
    LD_LIBRARY_PATH=$OPT/lib \
    LC_ALL=en_US.UTF-8 \
    LANG=en_US.UTF-8

WORKDIR /tmp/build

# build tools from other repos
COPY build/opt-build.sh build/
RUN bash build/opt-build.sh $OPT

COPY build/opt-build-local.sh build/
COPY c++ c++
COPY perl perl

# build the tools in this repo, separate to reduce build time on errors
RUN bash build/opt-build-local.sh $OPT

FROM ubuntu:20.04

LABEL maintainer="cgphelp@sanger.ac.uk" \
      uk.ac.sanger.cgp="Cancer, Ageing and Somatic Mutation, Wellcome Trust Sanger Institute" \
      description="cgpPindel docker"

# hadolint ignore=DL3008
RUN apt-get -yq update \
&& apt-get install -yq --no-install-recommends \
apt-transport-https \
locales \
curl \
ca-certificates \
libperlio-gzip-perl \
bzip2 \
psmisc \
time \
zlib1g \
liblzma5 \
libncurses5 \
p11-kit \
libcurl3-gnutls \
libcurl4 \
moreutils \
google-perftools \
unattended-upgrades && \
unattended-upgrade -d -v && \
apt-get remove -yq unattended-upgrades && \
apt-get autoremove -yq \
&& rm -rf /var/lib/apt/lists/* \
&& locale-gen en_US.UTF-8 \
&& update-locale LANG=en_US.UTF-8

ENV OPT /opt/wtsi-cgp
ENV PATH=$OPT/bin:$OPT/biobambam2/bin:$PATH \
    PERL5LIB=$OPT/lib/perl5 \
    LD_LIBRARY_PATH=$OPT/lib \
    LC_ALL=en_US.UTF-8 \
    LANG=en_US.UTF-8

COPY --from=builder $OPT $OPT

## USER CONFIGURATION
RUN adduser --disabled-password --gecos '' ubuntu && chsh -s /bin/bash && mkdir -p /home/ubuntu

USER    ubuntu
WORKDIR /home/ubuntu

CMD ["/bin/bash"]
