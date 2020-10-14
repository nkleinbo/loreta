FROM ubuntu:18.04

LABEL maintainer="nkleinbo@cebitec.uni-bielefeld.de"
LABEL version="0.1"
LABEL description="This is a custom Docker image to run \
the T-DNA Nanopore Analysis tool, see: https://github.com/nkleinbo/tdna_nanopore"

ARG DEBIAN_FRONTEND=noninteractive

RUN apt update

#Install python3, pip and Pillow:
RUN apt install -y python3 && \
    apt install -y python3-pip && \
    pip3 install Pillow

#Install gnu parallel 
RUN apt install -y parallel

#Install samtools
RUN apt install -y bzip2 && \
    apt install -y zlib1g-dev && \
    apt install -y libncurses5-dev && \
    apt install -y libbz2-dev && \
    apt install -y liblzma-dev && \
    apt install -y libcurl4-gnutls-dev && \
    apt install -y libssl-dev && \
    apt install -y build-essential && \
    apt install -y wget && \
    wget https://github.com/samtools/samtools/releases/download/1.11/samtools-1.11.tar.bz2 && \
    tar -xf samtools-1.11.tar.bz2 && \
    rm samtools-1.11.tar.bz2 && \
    cd samtools-1.11 && \
    ./configure && \
    make && \
    make install && \
    cd .. && \
    rm -rf samtools-1.11

#Install minimap2
RUN apt install -y libz-dev && \
    apt install -y build-essential && \
    apt install -y git && \
    git clone --branch v2.17 https://github.com/lh3/minimap2 && \
    cd minimap2 && \
    make && \
    mv minimap2 /usr/local/bin/. && \
    cd .. && \
    rm -rf minimap2

#Install seqtk
RUN apt install -y build-essential && \
    apt install -y git && \
    git clone --branch v1.3 https://github.com/lh3/seqtk && \
    cd seqtk && \
    make && \
    mv seqtk /usr/local/bin/. && \
    cd .. && \
    rm -rf seqtk

#Install canu
RUN apt install -y build-essential && \
    apt install -y git && \
    git clone --branch v2.1 https://github.com/marbl/canu && \
    cd canu/src && \
    make && \
    cd ../../ && \
    mv canu /usr/local/lib/ && \
    ln -s /usr/local/lib/canu/build/bin/canu /usr/local/bin/canu

#Install EMBOSS
RUN apt install -y emboss

#Install bedtools
RUN apt install -y python && \
    apt install -y build-essential && \
    apt install -y git && \
    git clone --branch v2.29.2 https://github.com/arq5x/bedtools2 && \
    cd bedtools2 && \
    make && \
    mv bin/* /usr/local/bin && \
    cd .. && \
    rm -rf bedtools2

#Install BLAST+
RUN wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.10.1/ncbi-blast-2.10.1+-x64-linux.tar.gz && \
    tar -xzf ncbi-blast-2.10.1+-x64-linux.tar.gz && \
    mv ncbi-blast-2.10.1+/bin/* /usr/local/bin && \
    rm -rf ncbi-blast-2.10.1+*

#Install tdna_nanopore
RUN apt install -y git && \
    git clone https://github.com/nkleinbo/tdna_nanopore/ && \
    gunzip tdna_nanopore/references/*.gz


