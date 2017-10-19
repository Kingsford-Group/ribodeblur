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
from deblur_utils import get_rlen_range_from_vblur, batch_build_ctrue, merge_profiles
from deblur_result_io import read_essentials, read_vblur, build_cobs_with_shifts
from footprint_hist_parser import parse_rlen_hist, get_cds_range, get_transcript_profiles

def merge_sparse_profiles(ctrue, cobs):
    """ merge sparse profiles in cobs into ctrue """
    cmix = {}
    for tid in cobs:
        if tid not in ctrue:
            tid_new = tid + " (no deblur)"
            cmix[tid_new] = cobs[tid]
        else:
            cmix[tid] = {}
            for rlen in cobs[tid]:
                if rlen in ctrue[tid]: 
                    cmix[tid][rlen] = ctrue[tid][rlen]
                else: 
                    cmix[tid][rlen] = cobs[tid][rlen]
    return cmix
                    
def construct_deblur_profiles(deblur_fname, vblur_fname, hist_fname, cds_range):
    """ construct deblurred profiles from deblur results """
    ptrue, eps = read_essentials(deblur_fname)
    b = read_vblur(vblur_fname)
    vrlen_min, vrlen_max = get_rlen_range_from_vblur(b)
    tlist = parse_rlen_hist(hist_fname)
    tprofile = get_transcript_profiles(tlist, cds_range, utr5_offset, utr3_offset)
    cobs_shift = build_cobs_with_shifts(tprofile, cds_range, utr5_offset, utr3_offset, vrlen_min, vrlen_max, klist)
    ctrue = batch_build_ctrue(ptrue, eps, cobs_shift)
    cmix = merge_sparse_profiles(ctrue, cobs_shift)
    ctrue_merge = { tid: merge_profiles(cmix[tid]) for tid in cmix }
    print("total transcripts with footprints in coding regions:", len(ctrue_merge))
    return ctrue_merge

def batch_build_Aprof(prof_dic, cds_range, utr5_offset, asite_offset):
    """ trim deblurred profiles to only include coding regions """
    aprof = {}
    for tid, prof in prof_dic.items():
        idx = tid.find(" (no deblur)")
        if idx != -1: 
            tid_name = tid[:idx]
        else:
            tid_name = tid
        cds_start, cds_end = cds_range[tid_name]
        istart = utr5_offset - asite_offset
        iend = istart + ( cds_end - cds_start )
        try: aprof[tid] = prof[istart: iend]
        except IndexError:
            print(tid, prof)
            exit(1)
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
        cds_range = get_cds_range(cds_fa)
        ctrue_merge = construct_deblur_profiles(eps_fname, vblur_fname, raw_hist, cds_range)
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
