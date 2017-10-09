# Ribodelur: a python package for estimating A-site locations of ribosomes from ribo-seq

## Overview
Ribodeblur is a suite of python scripts that perform noise reduction for ribosome profiling data. Ribosome profiling is a high-throughput-sequencing technique that captures mRNA fragments underneath ribosomes during translation. It is widely used to study translational dynamics and mechanics. Despite the fixed size of a ribosome, the size of captured mRNA fragments can span a wide range due to imperfect digestion during the complicated sequencing process. Such a side effect distort ribosome locations and thus can confound interpretations of the translational process. We developed a computational method to recover a clear ribiosome A-site position signal from ribosome profiling data. Our method successfully reduce noise from ribosome profiles and we hope our method can advance discoveries of translational regulations.

## Citation
This is the package for our previous work appeared at RECOMB in 2016: 

__Hao Wang, Joel McManus and Carl Kingsford__. *Accurate recovery of ribosome position signals reveals slow translation of wobble-pairing codons in yeast*. In RECOMB 2016 and JCB 2016.

## Prerequisites
`ribodeblur` is a python package based on `python 3`, and all of the prerequisites of `ribodeblur` are widely used python packages. They can be easily installed with [`easy_install`](http://peak.telecommunity.com/DevCenter/EasyInstall#using-easy-install),[`pip`](https://pip.pypa.io/en/stable/user_guide/) or [`conda`](https://conda.io/docs/user-guide/getting-started.html). 
* python (3.6)
* python packages: `numpy (1.13.0)`, `scipy (0.19.1)`, `bcbiogff (0.6.4)`, `biopython (1.68)`, and `pysam (0.11.2.2)`.
* STAR (2.5.3a) for aligning ribo-seq reads to the transcriptome.

## Usage
There are three major steps for running `ribodeblur`. First `ribodeblur` prepares a transcriptome `fasta` file that includes UTR regions for all included transcripts. Second, `ribodeblur` maps ribo-seq reads of a given sample to the transcriptome with STAR. Third, `ribodeblur` generates A-site profiles given the transcriptome reference and the read alignments. 

### Step 1: generate transcriptome fasta file
To generate a transcriptome fasta file, run:
```
python build_reference.py -f GENOME.FA -a GENOME.GFF -o TRANSCRIPTOME.FA
```
where `GENOME.FA` is the genome fasta file and `GENOME.GFF` is the annotation file, both should be downloaded from [Ensembl](https://www.ensembl.org/info/data/ftp/index.html). `TRANSCRIPTOME.FA` is the transcriptome reference the script generated. 

parameter `-p` specifiies the UTR padding regions added to each transcript, default is 100.

the output fasta file looks like this:
<blockquote>
>YAL069W |CDS:100-415|length:515
CATATCCAACCCACTGCCACTTACCCTACCATTACCCTACCATCCACCATGACCTACTCA
CCATACTGTTCTTCTACCCACCATATTGAAACGCTAACAAATGATCGTAAATAACACACA
...
</blockquote>

Here the CDS range provided in the header is important. `ribodeblur` relies on these ranges for down-stream analysis.

### Step 2: align ribo-seq reads to the transcriptome
To align ribo-seq reads to the transcriptome, run:
```
python map_to_reference.py -c CONTAMINANT.FA -t TRANCRITPOME.FA -r RIBOSEQ.FQ.GZ -id STAR_IDX_DIR -ad STAR_ALIGN_DIR
```
where `CONTAMINANT.FA` is the contaminant reference fasta (e.g. rRNA, tRNA, etc.), `TRANSCRIPTOME.FA` is the transcriptome reference generated from the previous step, `RIBOSEQ.FQ.GZ` is the ribo-seq raw reads file in either fasta or fastq format (`.gz` also supported), `STAR_IDX_DIR` is the directory to store STAR index, and `STAR_ALIGN_DIR` is the directory to store the alignment results (in `BAM` format).

The final output of this step is an alignment file ends with `_transcript_Aligned.out.bam` under directory `STAR_ALIGN_DIR`. 


