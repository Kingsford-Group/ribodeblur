#!/bin/bash

# download yeast reference genome and annotation from Ensembl
# inputs
work_dir="$HOME/data/ribodeblur/refs"
genome_link="ftp://ftp.ensembl.org/pub/release-89/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz"
gff_link="ftp://ftp.ensembl.org/pub/release-89/gff3/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.89.gff3.gz"
ncrna_link="ftp://ftp.ensembl.org/pub/release-89/fasta/saccharomyces_cerevisiae/ncrna/Saccharomyces_cerevisiae.R64-1-1.ncrna.fa.gz"

mkdir -p ${work_dir}
pushd ${work_dir}
curl -O ${genome_link}
curl -O ${gff_link}
curl -O ${ncrna_link}
gunzip -k -f *.gz
popd
