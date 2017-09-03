#!/bin/bash

# download yeast riboseq data from SRA
# inputs
work_dir="$HOME/data/ribodeblur/data"
riboseq_link="ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX476%2FSRX476345/SRR1177157/SRR1177157.sra"

mkdir -p ${work_dir}
pushd ${work_dir}
wget -P ${work_dir} -N ${riboseq_link}
riboseq_sra=`basename ${riboseq_link}`
fastq-dump ${work_dir}/${riboseq_sra} --gzip
popd
