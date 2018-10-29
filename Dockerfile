# Base image
FROM ubuntu:16.04

################## METADATA ######################

LABEL base_image="ubuntu:16.04"
LABEL version="4"
LABEL software="Biocontainers base Image"
LABEL software.version="1.0.0"
LABEL about.summary="Base image for BioDocker"
LABEL about.home="http://biocontainers.pro"
LABEL about.documentation="https://github.com/BioContainers/specs/wiki"
LABEL about.license_file="https://github.com/BioContainers/containers/blob/master/LICENSE"
LABEL about.license="SPDX:Apache-2.0"
LABEL about.tags="Genomics,Proteomics,Transcriptomics,General,Metabolomics"

################## MAINTAINER ######################
MAINTAINER Felipe da Veiga Leprevost <felipe@leprevost.com.br>

ENV DEBIAN_FRONTEND noninteractive

RUN mv /etc/apt/sources.list /etc/apt/sources.list.bkp && \
    bash -c 'echo -e "deb mirror://mirrors.ubuntu.com/mirrors.txt xenial main restricted universe multiverse\n\
deb mirror://mirrors.ubuntu.com/mirrors.txt xenial-updates main restricted universe multiverse\n\
deb mirror://mirrors.ubuntu.com/mirrors.txt xenial-backports main restricted universe multiverse\n\
deb mirror://mirrors.ubuntu.com/mirrors.txt xenial-security main restricted universe multiverse\n\n" > /etc/apt/sources.list' && \
    cat /etc/apt/sources.list.bkp >> /etc/apt/sources.list && \
    cat /etc/apt/sources.list

RUN apt-get clean all && \
    apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y  \
        autotools-dev   \
        automake        \
        cmake           \
        curl            \
        grep            \
        sed             \
        dpkg            \
        fuse            \
        git             \
        wget            \
        zip             \
        openjdk-8-jre   \
        build-essential \
        pkg-config      \
        python          \
	python-dev      \
        python-pip      \
        bzip2           \
        ca-certificates \
        libglib2.0-0    \
        libxext6        \
        libsm6          \
        libxrender1     \
        git             \
        mercurial       \
        subversion      \
        zlib1g-dev &&   \
        apt-get clean && \
        apt-get purge && \
        rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN echo 'export PATH=/opt/conda/bin:$PATH' > /etc/profile.d/conda.sh && \
    wget --quiet https://repo.continuum.io/miniconda/Miniconda2-4.0.5-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh

RUN TINI_VERSION=`curl https://github.com/krallin/tini/releases/latest | grep -o "/v.*\"" | sed 's:^..\(.*\).$:\1:'` && \
    curl -L "https://github.com/krallin/tini/releases/download/v${TINI_VERSION}/tini_${TINI_VERSION}.deb" > tini.deb && \
    dpkg -i tini.deb && \
    rm tini.deb && \
    apt-get clean

RUN mkdir /data /config

# give write permissions to conda folder
RUN chmod 777 -R /opt/conda/

ENV PATH=$PATH:/opt/conda/bin

RUN conda config --add channels r
RUN conda config --add channels bioconda

RUN conda upgrade conda

RUN conda install bwa=0.7.15

USER root
WORKDIR /

# Install SRA-Toolkit
RUN apt-get update -y && apt-get install -y wget && \
    apt-get clean && rm -rf /var/cache/apt/archives/* /var/lib/apt/lists/* && \
    wget -P /usr/bin "https://raw.githubusercontent.com/inutano/pfastq-dump/master/bin/pfastq-dump" && \
    chmod +x /usr/bin/pfastq-dump && \
    wget -P / "http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.9.2/sratoolkit.2.9.2-ubuntu64.tar.gz" && \
    tar zxf /sratoolkit.2.9.2-ubuntu64.tar.gz && \
    cp -r /sratoolkit.2.9.2-ubuntu64/bin/* /usr/bin && \
    rm -fr /sratoolkit.2.9.2-ubuntu64*

WORKDIR /opt

RUN apt-get update && apt-get install -y \
        wget \
        unzip \
        build-essential \
        apt-utils --yes --force-yes \
        openssh-client \
        curl \
        zlib1g-dev \
        libncurses5-dev \
        git

# Install aspc (aspera-client)
ADD http://download.asperasoft.com/download/sw/ascp-client/3.5.4/ascp-install-3.5.4.102989-linux-64.sh /opt/

# No https, so verify sha1
RUN test $(sha1sum /opt/ascp-install-3.5.4.102989-linux-64.sh |cut -f1 -d\ ) = a99a63a85fee418d16000a1a51cc70b489755957 && \
    sh /opt/ascp-install-3.5.4.102989-linux-64.sh


# Install Aspera connect so we can get the private key
ADD https://download.asperasoft.com/download/sw/connect/3.8.1/ibm-aspera-connect-3.8.1.161274-linux-g2.12-64.tar.gz /opt/

WORKDIR /opt/

RUN mkdir -p /home/data/.aspera/connect
RUN chmod -R 777 /home/data/.aspera/connect/
RUN tar -xvf ibm-aspera-connect-3.8.1.161274-linux-g2.12-64.tar.gz
RUN useradd data
USER data

RUN ./ibm-aspera-connect-3.8.1.161274-linux-g2.12-64.sh

# Must switch back to root as user in order to run vdb-config properly
USER root

### Now private key can be found at: /home/data/.aspera/connect/etc/asperaweb_id_dsa.openssh

WORKDIR /

# set the environment variables
ENV samtools_version 1.9
ENV bcftools_version 1.9
ENV htslib_version 1.9

# run update and install necessary packages
RUN apt-get update -y && apt-get install -y \
    bzip2 \
    build-essential \
    zlib1g-dev \
    libncurses5-dev \
    libncursesw5-dev \
    libnss-sss \
    libbz2-dev \
    liblzma-dev

# download the suite of tools
ADD https://github.com/samtools/samtools/releases/download/${samtools_version}/samtools-${samtools_version}.tar.bz2 /usr/bin/
ADD https://github.com/samtools/bcftools/releases/download/${bcftools_version}/bcftools-${bcftools_version}.tar.bz2 /usr/bin/
ADD https://github.com/samtools/htslib/releases/download/${htslib_version}/htslib-${htslib_version}.tar.bz2 /usr/bin/

# extract files for the suite of tools
RUN tar -xjf /usr/bin/samtools-${samtools_version}.tar.bz2 -C /usr/bin/
RUN tar -xjf /usr/bin/bcftools-${bcftools_version}.tar.bz2 -C /usr/bin/
RUN tar -xjf /usr/bin/htslib-${htslib_version}.tar.bz2 -C /usr/bin/

# run make on the source
RUN cd /usr/bin/htslib-${htslib_version}/ && ./configure
RUN cd /usr/bin/htslib-${htslib_version}/ && make
RUN cd /usr/bin/htslib-${htslib_version}/ && make install

RUN cd /usr/bin/samtools-${samtools_version}/ && ./configure
RUN cd /usr/bin/samtools-${samtools_version}/ && make
RUN cd /usr/bin/samtools-${samtools_version}/ && make install

RUN cd /usr/bin/bcftools-${bcftools_version}/ && make
RUN cd /usr/bin/bcftools-${bcftools_version}/ && make install