# Ribodeblur: a python package for estimating A-site locations of ribosomes from ribo-seq

## Overview
Ribodeblur is a suite of python scripts that perform noise reduction for ribosome profiling data. Ribosome profiling is a high-throughput-sequencing technique that captures mRNA fragments underneath ribosomes during translation. It is widely used to study translational dynamics and mechanics. Despite the fixed size of a ribosome, the size of captured mRNA fragments can span a wide range due to imperfect digestion during the complicated sequencing process. Such a side effect distort ribosome locations and thus can confound interpretations of the translational process. We developed a computational method to recover a clear ribiosome A-site position signal from ribosome profiling data. Our method successfully reduce noise from ribosome profiles and we hope our method can advance discoveries of translational regulations.

## Citation
This is the package for our previous work appeared at RECOMB in 2016: 

__Hao Wang, Joel McManus and Carl Kingsford__. *Accurate recovery of ribosome position signals reveals slow translation of wobble-pairing codons in yeast*. In RECOMB 2016 and JCB 2016.

## Prerequisites
`ribodeblur` is a python package based on `python3`, and all prerequisites (except for STAR) are widely used python packages. They can be easily installed with [`pip`](https://pip.pypa.io/en/stable/user_guide/) or [`conda`](https://conda.io/docs/user-guide/getting-started.html). 
* python (3.6)
* python packages: `numpy (1.13.0)`, `scipy (0.19.1)`, `bcbiogff (0.6.4)`, `biopython (1.68)`, and `pysam (0.11.2.2)`.
* STAR (2.5.3a) for aligning ribo-seq reads to the transcriptome.

## Set up prerequisites with conda
[`Conda`](https://conda.io/docs/) is an environment and package management system. If `conda` is not installed in your system, follow [these](https://conda.io/docs/user-guide/install/index.html) instructions. Here lists key commands to setup `conda` on a linux machine:
* download and install `miniconda`:
```
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh
bash ~/miniconda.sh -b -p $HOME/miniconda
export PATH="$HOME/miniconda/bin:$PATH"
```
* create virtual environment with all prerequisites installed:
```
conda env create -f ribodeblur.yml
```
Here [`ribodeblur.yml`](https://github.com/Kingsford-Group/ribodeblur/blob/master/ribodeblur.yml) is the conda enviroment yaml file provided with the package. 
* activate conda environment:
```
source activate ribodeblur
```

## Usage
There are three major steps for running `ribodeblur`. First `ribodeblur` prepares a transcriptome `fasta` file that includes UTR regions for all transcripts. Second, `ribodeblur` maps ribo-seq reads of a given sample to the transcriptome with STAR. Third, `ribodeblur` generates A-site profiles given the transcriptome reference and the read alignments. 

### Step 1: generate transcriptome fasta file
The core of `ribodeblur` depends on reconcentrating off-frame ribo-seq reads back to their in-frame locations based on a blind-deconvolution approach. To achieve that, `ribodeblur` requires ribosome profiles from both coding regions and a small portion of untranslated regions.   

To generate a transcriptome fasta file with UTR regions included, run:
```
python build_reference.py -f GENOME.FA -a GENOME.GFF -o TRANSCRIPTOME.FA
```
where `GENOME.FA` is the genome fasta file and `GENOME.GFF` is the annotation file, both should be downloaded from [Ensembl](https://www.ensembl.org/info/data/ftp/index.html). `TRANSCRIPTOME.FA` is the transcriptome reference the script generated. 

Optional parameter `-p` specifies the UTR padding regions added to each transcript, default is 100.

the output fasta file looks like this:
<blockquote>
>YAL069W |CDS:100-415|length:515
CATATCCAACCCACTGCCACTTACCCTACCATTACCCTACCATCCACCATGACCTACTCA
CCATACTGTTCTTCTACCCACCATATTGAAACGCTAACAAATGATCGTAAATAACACACA
...
</blockquote>

Here the CDS range provided in the header is important. `ribodeblur` relies on these ranges for down-stream analysis.

### Step 2: align ribo-seq reads to the transcriptome
`ribodeblur` then aligns the raw ribo-seq reads to the transcriptome with STAR. It automatically builds STAR index of the reference, filters potential non-mRNA contaminants, and aligns the remaining reads to the transcriptome.

To align ribo-seq reads to the transcriptome, run:
```
python map_to_reference.py -c CONTAMINANT.FA -t TRANCRITPOME.FA -r RIBOSEQ.FQ -id STAR_IDX_DIR -ad STAR_ALIGN_DIR
```
where `CONTAMINANT.FA` is the contaminant reference fasta (e.g. rRNA, tRNA, etc.), `TRANSCRIPTOME.FA` is the transcriptome reference generated from **Step 1**, `RIBOSEQ.FQ` is the ribo-seq raw reads (in fasta/fastq format, also support `.gz`), `STAR_IDX_DIR` is the directory to store STAR index, and `STAR_ALIGN_DIR` is the directory to store the alignment results (in `BAM` format).

##### Optional parameters:
* `-p/--nproc`: number of threads to run STAR (default=30)
* `a/--adapter`: adapter sequence to trim off by STAR (default=`CTGTAGGCACCATCAAT`)
* `-f/--force`: force re-run of STAR alignment.

The final output of this step is an alignment file ends with `_transcript_Aligned.out.bam` under directory `STAR_ALIGN_DIR`. 

### Step 3: deblur ribosome profiles and output A-site profiles
The final step of `ribodeblur` is the _deblur_ step, where the read alignments are processed and grouped by read length for each transcript, and only profiles with high coverage (>50% non-zero loci) are kept. Then for each read length, a _blur vector_ is learned from a meta-profile combining all transcripts, and each transcript profile are _deblurred_ independently. Lastly, the deblurred profiles from different read length are merged together to get a final ribosome A-site profile.

To generate A-site profiles, run:
```
python deblur_pipeline.py -r TRANSCRIPTOME.FA -b ALIGN.BAM -o OUTPUT_PREFIX
```
where `TRANSCRITPOME.FA` is the transcriptome reference generated from **Step 1**, `ALIGN.BAM` is the alignment file generated form **Step 2**, and `OUTPUT_PREFIX` is the prefix of the output for the A-site profile.

The final otuput of this step is an A-site profile starts with `OUTPUT_PREFIX` and ends with `.profile`.
It looks like this:
<blockquote>
YBR248C: 3 0 0 1 0 0 0 0 0 1 0 0 132 0 15 186 1 0 171 ...
</blockquote>
Each line starts with the transcript name followed by the read count for each nucleotide location.

## Test case
`ribodeblur` requires the tested ribo-seq data with deep coverage. Here provides a real-world example from Albert _et al._

To download the yeast genome and non-coding reference (as contaminant reference), run script `download_refs.sh`. This automatically downloads the reference sequences from Ensembl to directory `$HOME/data/ribodeblur/refs`.

To download the ribo-seq data, run script `download_ribobseq.sh`. This downloads the test-case ribo-seq data to directory `$HOME/data/ribodeblur/data`.

To run the `ribodeblur` pipeline on the test set, use script `run_ribodeblur.sh`. This runs `ribodeblur` on the test set and outputs the A-site profile to `$Home/data/ribodeblur/deblur/by_fp.profile`.