# Ribodelur: a python package for estimating A-site locations of ribosomes from ribo-seq

## Overview
Ribodeblur is a suite of python scripts that perform noise reduction for ribosome profiling data. Ribosome profiling is a high-throughput-sequencing technique that captures mRNA fragments underneath ribosomes during translation. It is widely used to study translational dynamics and mechanics. Despite the fixed size of a ribosome, the size of captured mRNA fragments can span a wide range due to imperfect digestion during the complicated sequencing process. Such a side effect distort ribosome locations and thus can confound interpretations of the translational process. We developed a computational method to recover a clear ribiosome A-site position signal from ribosome profiling data. Our method successfully reduce noise from ribosome profiles and we hope our method can advance discoveries of translational regulations.

## Citation
This is the package for our previous work appeared at RECOMB in 2016: 

__Hao Wang, Joel McManus and Carl Kingsford__. *Accurate recovery of ribosome position signals reveals slow translation of wobble-pairing codons in yeast*. In RECOMB 2016 and JCB 2016.

## Prerequisites
`ribodeblur` is a python package based on `python 3`, and all of the prerequisites of `ribodeblur` are widely used python packages. They can be easily installed with [`easyinstall`](http://peak.telecommunity.com/DevCenter/EasyInstall#using-easy-install),[`pip`](https://pip.pypa.io/en/stable/user_guide/) or [`conda`](https://conda.io/docs/user-guide/getting-started.html). 
* python (3.6)
* python packages: `numpy (1.13.0)`, `scipy (0.19.1)`, `bcbiogff (0.6.4)`, `biopython (1.68)`, and `pysam (0.11.2.2)`.