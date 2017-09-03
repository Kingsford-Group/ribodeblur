#!/usr/bin/env python
"""
map reads to reference with STAR
step 0: build star index
step 1: filter contaminants
step 2: align to transcriptome
"""

import os
import sys
import argparse
from string import Template
import subprocess

import numpy as np
import Bio.SeqIO

contaminant_idx_cmd = Template("""STAR --runMode genomeGenerate --genomeFastaFiles ${ref_fa} --genomeDir ${ref_idx} --runThreadN ${nproc} ${param_str}""")
common_align_param = Template("""--runThreadN ${nproc} --clip3pAdapterSeq ${adapter} --seedSearchLmax 10 --outFilterMultimapScoreRange 0 --outFilterMultimapNmax 255 --outFilterMismatchNmax 1 --outFilterIntronMotifs RemoveNoncanonical """)
algin_con_cmd = Template("""STAR --genomeDir ${idx_dir} --readFilesIn ${in_fq} --outFileNamePrefix ${out_prefix} --outStd SAM --outReadsUnmapped Fastx --outSAMmode NoQS ${common_param} > /dev/null""")
align_ref_cmd = Template("""STAR --genomeDir ${idx_dir} --readFilesIn ${in_fq} --outFileNamePrefix ${out_prefix} --outSAMtype BAM Unsorted --outSAMmode NoQS --outSAMattributes NH NM ${common_param}""")

idx_files = ( "chrLength.txt", "chrNameLength.txt", "chrName.txt",
              "chrStart.txt", "Genome", "genomeParameters.txt",
              "SA", "SAindex" )

def get_file_core(file_name):
    """ strip off folder and extension of a file name """
    basename = os.path.basename(file_name)
    noext1 = os.path.splitext(basename)[0]
    noext2 = os.path.splitext(noext1)
    if noext2[1] in set([".fq", ".fa", ".fasta", ".fastq", ".FQ", ".FA", ".FASTA", ".FASTQ"]):
        return noext2[0]
    else:
        return noext1

def get_ref_len_cnt(fa_fname):
    """ get total reference length and counts for fasta file """
    ref_cnt = 0
    ref_len = 0
    for record in Bio.SeqIO.parse(fa_fname, "fasta"):
        ref_cnt += 1
        ref_len += len(record.seq)
    return ref_len, ref_cnt

def set_ref_params(ref_len, ref_cnt):
    """ 
    set star index params as described in star manual 
    https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf
    page 6
    section 2.2.5 Very small genome
    section 2.2.6 Genome with a large number of references
    """
    # --genomeSAindexNbases
    SAI_bits = min(14, np.ceil(np.log2(ref_len)/2 - 1) )
    cmd = "--genomeSAindexNbases {:.0f} ".format(SAI_bits)
    # --genomeChrBinNbits
    if ref_cnt > 5000:
        chr_bits = min(18, np.ceil(np.log2(ref_len/ref_cnt)))
        cmd += "--genomeChrBinNbits {:.0f} ".format(chr_bits)
    return cmd

def check_files(dir_name, file_list):
    """ check whether all listed files exists in a folder """
    for fname in file_list:
        if not os.path.exists("{}/{}".format(dir_name, fname)):
            return False
    return True

def build_star_idx(ref_fa, star_idx_dir, nproc, force=True):
    """ build STAR index of a given fasta """
    # check inputs
    if not os.path.exists(ref_fa): 
        print("FATAL: reference fasta {} not exist!".format(ref_fa), file=sys.stderr)
        print("abort program!", file=sys.stderr)
        exit(1)
    # prepare inputs
    idx_dir = "{}/{}".format(star_idx_dir, get_file_core(ref_fa))
    if not os.path.exists(idx_dir): os.makedirs(idx_dir)
    if not check_files(idx_dir, idx_files) or force:
        # prepare parameters
        ref_len, ref_cnt = get_ref_len_cnt(ref_fa)
        param_str = set_ref_params(ref_len, ref_cnt)
        cmd = contaminant_idx_cmd.substitute(ref_fa = ref_fa,
                                             ref_idx = idx_dir,
                                             nproc = 30,
                                             param_str = param_str)
        print(cmd, file=sys.stderr)
        # run command
        subprocess.call(cmd, shell=True)
        # after-run check
        if not check_files(idx_dir, idx_files):
            print("FATAL: map_to_reference failed at building index for {}!".format(ref_fa), file=sys.stderr)
            print("abort program!", file=sys.stderr)
            exit(1)
    print("Done building index for {}, output to {}.".format(ref_fa, idx_dir), file=sys.stderr)
    return idx_dir

def filter_contaminant(star_idx, read_fq, adapter, nproc, align_dir, force=True):
    """ run STAR to filter contaminants """
    # check inputs
    if not os.path.exists(read_fq): 
        print("FATAL: read file {} not exist!".format(read_fq), file=sys.stderr)
        print("abort program!", file=sys.stderr)
        exit(1)
    # prepare inputs
    if not os.path.exists(align_dir): os.makedirs(align_dir)
    out_prefix = "{}/{}_rrna_".format(align_dir, get_file_core(read_fq))
    filter_fq = "{}Unmapped.out.mate1".format(out_prefix)
    if not os.path.exists(filter_fq) or force:
        # prepare parameters
        common_param = common_align_param.substitute(nproc = nproc, adapter = adapter)
        if read_fq.endswith(".gz"): common_param += '--readFilesCommand "zcat <"'
        cmd = algin_con_cmd.substitute(idx_dir = star_idx,
                                       in_fq = read_fq,
                                       out_prefix = out_prefix,
                                       common_param = common_param)
        print(cmd, file=sys.stderr)
        # run command
        subprocess.call(cmd, shell=True)
        # after-run check
        if not os.path.exists(filter_fq):
            print("FATAL: map_to_reference failed at filtering reads {}!".format(read_fq), file=sys.stderr)
            print("abort program!", file=sys.stderr)
            exit(1)
    print("Done filtering reads {}, output to {}.".format(read_fq, filter_fq), file=sys.stderr)
    return filter_fq

def align_to_transcript(star_idx, filter_fq, adapter, nproc, force=True):
    """ run STAR to align remaining reads to the transcriptome """
    # check inputs
    if not os.path.exists(filter_fq): 
        print("FATAL: filtered reads {} not exist!".format(filter_fq), file=sys.stderr)
        print("abort program!", file=sys.stderr)
        exit(1)
    # prepare inputs
    out_prefix = filter_fq.rstrip("_rrna_Unmapped.out.mate1")+"_transcript_"
    map_bam = out_prefix + "Aligned.out.bam"
    if not os.path.exists(map_bam) or force:
        # prepare parameters
        common_param = common_align_param.substitute(nproc = nproc, adapter = adapter)
        cmd = align_ref_cmd.substitute(idx_dir = star_idx,
                                       in_fq = filter_fq,
                                       out_prefix = out_prefix,
                                       common_param = common_param)
        print(cmd, file=sys.stderr)
        # run command
        subprocess.call(cmd, shell=True)
        # after-run check
        if not os.path.exists(map_bam):
            print("FATAL: map_to_reference failed at mapping reads {}!".format(filter_fq), file=sys.stderr)
            print("abort program!", file=sys.stderr)
            exit(1)
    print("Done mapping to transcriptome for {}, output to {}.".format(filter_fq, map_bam), file=sys.stderr)
    return map_bam
    
def make_arg_parser():
    """ command line arguments """
    parser = argparse.ArgumentParser(prog="map_to_reference.py", add_help=True, description="map ribo-seq reads to the transcriptome with STAR")
    parser.add_argument("-c", "--contaminant", required=True, help="reference fasta of potential contaminants (e.g. rRNA, tRNA, etc.)")
    parser.add_argument("-t", "--transcript", required=True, help="reference fasta of padded transcriptome")
    parser.add_argument("-r", "--riboseq", required = True, help="raw ribo-seq reads (in fasta/fastq, can be compressed -- has to end with '.gz')")
    parser.add_argument("-id", "--idx_dir", required = True, help="directory for storing STAR index")
    parser.add_argument("-ad", "--align_dir", required = True, help="directory for storing STAR alignments")
    parser.add_argument("-p", "--nproc", type=int, default=30, help="number of threads to run STAR (default=30)")
    parser.add_argument("-a", "--adapter", default="CTGTAGGCACCATCAAT", help="adapter sequence to trim off by STAR (default='CTGTAGGCACCATCAAT')")
    parser.add_argument("-f", "--force", action="store_true", default=False, help="force re-run of entire pipeline if provided")
    return parser

def main():
    """ pipeline script to call STAR for alignment ribo-seq reads to the transcriptome """
    parser = make_arg_parser()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    print("building contaminant index...", file=sys.stderr)
    con_idx_dir = build_star_idx(args.contaminant, args.idx_dir, args.nproc, args.force)
    print("building transcriptome index...", file=sys.stderr)
    ref_idx_dir = build_star_idx(args.transcript, args.idx_dir, args.nproc, args.force)
    print("filterint contaminants...", file=sys.stderr)
    filter_fq = filter_contaminant(con_idx_dir, args.riboseq, args.adapter, args.nproc, args.align_dir, args.force)
    print("mapping to transcriptome...", file=sys.stderr)
    map_bam = align_to_transcript(ref_idx_dir, filter_fq, args.adapter, args.nproc, args.force)

if __name__ == "__main__": main()
# python map_to_reference.py -c ~/data/ribodeblur/refs/Saccharomyces_cerevisiae.R64-1-1.ncrna.fa -t ~/data/ribodeblur/refs/Saccharomyces_cerevisiae.R64-1-1.transcript_100.fa -r ~/data/ribodeblur/data/SRR1177157.fastq.gz -id ~/data/ribodeblur/star_index -ad ~/data/ribodeblur/star_align
