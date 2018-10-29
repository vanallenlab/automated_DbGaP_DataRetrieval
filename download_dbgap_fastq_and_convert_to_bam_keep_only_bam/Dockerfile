################## BASE IMAGE ######################
FROM biocontainers/biocontainers:latest

################## METADATA ######################
LABEL base_image="biocontainers:latest"
LABEL version="2"
LABEL software="bwa"
LABEL software.version="0.7.15"
LABEL about.summary="Burrow-Wheeler Aligner for pairwise alignment between DNA sequences"
LABEL about.home="http://bio-bwa.sourceforge.net/"
LABEL about.documentation="http://bio-bwa.sourceforge.net/"
LABEL about.license_file="http://bio-bwa.sourceforge.net/"
LABEL about.license="SPDX:GPL-3.0"
LABEL extra.identifiers.biotools="bwa"
LABEL about.tags="Genomics"

################## MAINTAINER ######################
MAINTAINER Saulo Alves Aflitos <sauloal@gmail.com>

################## INSTALLATION ######################

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