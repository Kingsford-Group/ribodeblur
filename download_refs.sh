#!/bin/bash

# download yeast reference genome and annotation from Ensembl
# inputs
work_dir="$HOME/data/ribodeblur/refs"
genome_link="ftp://ftp.ensembl.org/pub/release-89/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz"
gff_link="ftp://ftp.ensembl.org/pub/release-89/gff3/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.89.gff3.gz"
ncrna_link="ftp://ftp.ensembl.org/pub/release-89/fasta/saccharomyces_cerevisiae/ncrna/Saccharomyces_cerevisiae.R64-1-1.ncrna.fa.gz"

# for sanity checks
cdna_link="ftp://ftp.ensembl.org/pub/release-89/fasta/saccharomyces_cerevisiae/cdna/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz"
cds_link="ftp://ftp.ensembl.org/pub/release-89/fasta/saccharomyces_cerevisiae/cds/Saccharomyces_cerevisiae.R64-1-1.cds.all.fa.gz"
pep_link="ftp://ftp.ensembl.org/pub/release-89/fasta/saccharomyces_cerevisiae/pep/Saccharomyces_cerevisiae.R64-1-1.pep.all.fa.gz"

mkdir -p ${work_dir}
pushd ${work_dir}
wget -P ${work_dir} -N ${genome_link}
wget -P ${work_dir} -N ${gff_link}
wget -P ${work_dir} -N ${ncrna_link}
wget -P ${work_dir} -N ${cdna_link}
wget -P ${work_dir} -N ${cds_link}
wget -P ${work_dir} -N ${pep_link}
gunzip -k -f *.gz
popd
