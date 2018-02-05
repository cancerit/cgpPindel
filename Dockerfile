# This Dockerfile imports from cancerit/cgp-vcf:2.2.1 wich already includes
# cgpBigWig, PCAP-core and cgpVcf. cgpPindel system dependencies are satisfied
# with a successful install of cgpVcf.

# The cancerit/cgp-vcf container includes the environment variable OPT,
# this is reused here to install cgpPindel.
# As such there is no need to update PATH and PERL5LIB.

# Locale is also set to:
# ENV LC_ALL en_US.UTF-8
# ENV LANG en_US.UTF-8
FROM cancerit/cgp-vcf:2.2.1

# Set maintainer labels.
LABEL maintainer Keiran M. Raine <kr2@sanger.ac.uk>

# Add repo.
COPY . /code

# Install package.
RUN \
    cd /code && \
    ./setup.sh $OPT && \
    cd ~ && \
    rm -rf /code

# Set volume to data as per:
# https://github.com/BD2KGenomics/cgl-docker-lib
VOLUME /data
WORKDIR /data
