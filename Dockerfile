FROM debian:8.5

MAINTAINER Migun Shakya <migun@lanl.gov>

LABEL version="1.0"
LABEL software="PhaME"
LABEL description="Phylogenetic Analysis and Molecular Evolution"
LABEL website="https://github.com/mshakya/PhaME-1"
LABEL tags="genomics"

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

RUN apt-get update --fix-missing && apt-get install -y wget bzip2 ca-certificates \
    libglib2.0-0 libxext6 libsm6 libxrender1 \
    git curl make gcc gfortran g++ grep sed

RUN mkdir -p /opt && \
	cd /opt && git clone https://github.com/mshakya/PhaME-1.git && \
    cd PhaME-1 && bash INSTALL.sh &&
    echo 'export PATH=/opt/PhaME-1:/opt/PhaME-1/ext/miniconda/bin:$PATH' > /etc/profile.d/PhaME-1.sh && \
    echo 'export R_LIBS=/opt/targetngs/ext/lib/R_libs' >> /etc/profile.d/PhaME-1.sh && \
    apt-get clean

ENV PATH PhaME-1/thirdParty/miniconda/bin/:$PATH

