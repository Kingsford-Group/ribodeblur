#!/usr/bin/env python
"""
pipeline to run deblur process given reference and read alignment
step 1: generate length-specific profiles
step 2: filter high-coverage profiles
step 3: train blur vector from meta profile
step 4: deblur high-coverage profiles
step 5: combine length-specific profiles
"""

import os
import sys
import argparse
import operator
from global_params import *
from group_reads_by_length import group_reads_by_length
from high_cover_profile import filter_high_cover_profile
from train_vblur_from_meta import train_vblur_from_meta
from deblur_transcripts import deblur_transcripts
from deblur_utils import batch_build_ctrue, merge_profiles
from deblur_result_io import read_essentials
from footprint_hist_parser import parse_rlen_hist, get_cds_range

def get_tot_cnt(prof):
    """ count total reads (input: [(pos, cnt)] """
    return sum(map(operator.itemgetter(1), prof))
    
def get_rlen_tot_cnt(tprof):
    """ count total reads per length (input: {rlen: prof}) """
    return { rlen: get_tot_cnt(tprof[rlen]) for rlen in tprof }

def get_tot_cnt_from_tlist(tlist):
    """ convert raw rlen hist parser to total read count """
    return { t["tid"] : get_rlen_tot_cnt(t["prof"]) for rid,t in tlist.items() }

def count_reads_per_len(hist_fname):
    """ count total number of reads per read length per transcript """
    tlist = parse_rlen_hist(hist_fname)
    return get_tot_cnt_from_tlist(tlist)
    
def construct_deblur_profiles(deblur_fname, hist_fname):
    """ construct deblurred profiles from deblur results """
    ptrue, eps = read_essentials(deblur_fname)
    tot_cnts = count_reads_per_len(hist_fname)
    ctrue = batch_build_ctrue(ptrue, eps, tot_cnts)
    ctrue_merge = { tid: merge_profiles(ctrue[tid]) for tid in ctrue }
    return ctrue_merge

def batch_build_Aprof(prof_dic, cds_range, utr5_offset, asite_offset):
    """ trim deblurred profiles to only include coding regions """
    aprof = {}
    for tid, prof in prof_dic.items():
        cds_start, cds_end = cds_range[tid]
        istart = utr5_offset - asite_offset
        iend = istart + ( cds_end - cds_start )
        aprof[tid] = prof[istart: iend]
    return aprof

def write_profiles(profiles, profile_fname):
    """ write merged profiles to file """
    text = []
    for tid in sorted(profiles.keys()):
        prof_str = " ".join([ "{:.0f}".format(profiles[tid][i]) 
                              for i in range(len(profiles[tid])) ])
        text.append("{}: {}".format(tid, prof_str))
    ofile = open(profile_fname, 'w')
    ofile.write("\n".join(text))
    ofile.close()

def deblur_pipeline(bam_fname, cds_fa, oprefix, force):
    """ full pipeline for deblur ribo profiles for a given sample """
    # step 0: prepare input parameters
    odir = os.path.dirname(oprefix)
    if odir and not os.path.exists(odir): os.makedirs(odir)
    if oprefix.endswith("/"): oprefix += "ribo"
    raw_hist = "{}_raw.hist".format(oprefix)
    high_cov_hist = "{}_hc.hist".format(oprefix)
    vblur_fname = "{}.vblur".format(oprefix)
    eps_fname = "{}.eps".format(oprefix)
    profile_fname = "{}.profile".format(oprefix)
    # step 1: generate length-specific profiles
    if not os.path.exists(raw_hist) or force == True:
        group_reads_by_length(bam_fname, raw_hist)
        if not os.path.exists(raw_hist) or os.path.getsize(raw_hist) == 0:
            print("FATAL: deblur_pipeline failed at group_reads_by_len!", file=sys.stderr)
            print("abort program!", file=sys.stderr)
            exit(1)
    else: print("length-specific profile exists, use cached", file=sys.stderr)
    # step 2: filter high-coverage profiles
    if not os.path.exists(high_cov_hist) or force == True:
        filter_high_cover_profile(raw_hist, cds_fa, cover_ratio, cnt_threshold, high_cov_hist)
        if not os.path.exists(high_cov_hist) or os.path.getsize(high_cov_hist) == 0:
            print("FATAL: deblur_pipeline failed at filter_high_cover_profile!", file=sys.stderr)
            print("abort program!", file=sys.stderr)
            exit(1)
    else: print("high-coverage profile exists, use cached", file=sys.stderr)
    # step 3: train blur vector from meta profiles
    if not os.path.exists(vblur_fname) or force == True:
        train_vblur_from_meta(high_cov_hist, cds_fa, vblur_fname)
        if not os.path.exists(vblur_fname) or os.path.getsize(vblur_fname) == 0:
            print("FATAL: deblur_pipeline failed at train_vblur_from_meta!", file=sys.stderr)
            print("abort program!", file=sys.stderr)
            exit(1)
    else: print("vblur file exists, use cached", file=sys.stderr)
    # step 4: deblur high-coverage profiles
    if not os.path.exists(eps_fname) or force == True:
        deblur_transcripts(high_cov_hist, cds_fa, vblur_fname, eps_fname)
        if not os.path.exists(eps_fname) or os.path.getsize(eps_fname) == 0:
            print("FATAL: deblur_pipeline failed at deblur_transcripts!", file=sys.stderr)
            print("abort program!", file=sys.stderr)
            exit(1)
    else: print("deblur file exits, use cached", file=sys.stderr)
    # step 5: combine length-specific profiles
    if not os.path.exists(profile_fname) or force == True:
        ctrue_merge = construct_deblur_profiles(eps_fname, high_cov_hist)
        cds_range = get_cds_range(cds_fa)
        aprof = batch_build_Aprof(ctrue_merge, cds_range, -utr5_offset, asite_offset) 
        write_profiles(aprof, profile_fname)
        if not os.path.exists(profile_fname) or os.path.getsize(profile_fname) == 0:
            print("FATAL: deblur_pipeline failed at combine_profile!", file=sys.stderr)
            print("abort program!", file=sys.stderr)
            exit(1)
    else: print("final results exists, nothing needs to be done", file=sys.stderr)

def make_arg_parser():
    """ command line arguments """
    parser = argparse.ArgumentParser(prog="deblur_pipeline.py", add_help=True, description="pipeline to run deblur process given reference and read alignment")
    parser.add_argument("-r", "--cds_fa", required=True, help="Reference fasta of padded transcriptome")
    parser.add_argument("-b", "--bam", required=True, help="Alignment BAM file that maps reads to transcriptome")
    parser.add_argument("-o", "--oprefix", required = True, help="Prefix of output files")
    parser.add_argument("-f", "--force", action="store_true", default=False, help="force re-run of entire pipeline if provided")
    return parser

def main():
    """ command line wrapper """
    parser = make_arg_parser()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    deblur_pipeline(args.bam, args.cds_fa, args.oprefix, args.force)

if __name__ == "__main__": main()
# python deblur_pipeline.py -r /home/hw1/data/ribodeblur/refs/Saccharomyces_cerevisiae.R64-1-1.transcript_100.fa -b /home/hw1/data/ribodeblur/star_align/SRR1177157_transcript_Aligned.out.bam -o /home/hw1/data/ribodeblur/deblur/by_fp    
