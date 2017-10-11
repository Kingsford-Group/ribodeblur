#!/bin/bash

# test case for running ribodeblur

# data name and locations
genome_fa="$HOME/data/ribodeblur/refs/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa"
genome_gff="$HOME/data/ribodeblur/refs/Saccharomyces_cerevisiae.R64-1-1.89.gff3"
contaminant_fa="$HOME/data/ribodeblur/refs/Saccharomyces_cerevisiae.R64-1-1.ncrna.fa"
ribo_fq="$HOME/data/ribodeblur/data/SRR1177157.fastq.gz"
cds_fa="$HOME/data/ribodeblur/refs/Saccharomyces_cerevisiae.R64-1-1.transcript_100.fa"
star_idx="$HOME/data/ribodeblur/star_index"
align_dir="$HOME/data/ribodeblur/star_align"
align_bam="${align_dir}/SRR1177157_transcript_Aligned.out.bam"
oprefix="$HOME/data/ribodeblur/deblur/by_fp"

# step 1: build reference
python build_reference.py -f ${genome_fa} -a ${genome_gff} -o ${cds_fa}

# step 2: map reads to reference
python map_to_reference.py -c ${contaminant_fa} -t ${cds_fa} -r ${ribo_fq} -id ${star_idx} -ad ${align_dir}

# step 3: deblur profiles
python deblur_pipeline.py -r ${cds_fa} -b ${align_bam} -o ${oprefix}
