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
import operator
from group_reads_by_length import group_reads_by_length
from high_cover_profile import filter_high_cover_profile
from train_vblur_from_meta import train_vblur_from_meta
from deblur_transcripts import deblur_transcripts
from deblur_utils import batch_build_ctrue, merge_profiles
from deblur_result_io import read_essentials
from footprint_hist_parser import parse_rlen_hist

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

if __name__ == "__main__":
    force = False
    bam_fname = "/home/hw1/data/ribodeblur/star_align/SRR1177157_transcript_Aligned.out.bam"
    cds_fa = "/home/hw1/data/ribodeblur/refs/Saccharomyces_cerevisiae.R64-1-1.transcript_100.fa"
    odir = "/home/hw1/data/ribodeblur/deblur"
    file_core = "by_fp"
    raw_hist = "{}/{}_raw.hist".format(odir, file_core)
    high_cov_hist = "{}/{}_hc.hist".format(odir, file_core)
    vblur_fname = "{}/{}.vblur".format(odir, file_core)
    eps_fname = "{}/{}.eps".format(odir, file_core)
    profile_fname = "{}/{}.profile".format(odir, file_core)
    cover_ratio = 0.5
    cnt_threshold = 0
    
    # step 1: generate length-specific profiles
    if not os.path.exists(raw_hist) or force == True:
        group_reads_by_len_pipe(bam_fname, raw_hist)
    else:
        print("length-specific profile exists, use cached", file=sys.stderr)
    # step 2: filter high-coverage profiles
    if not os.path.exists(high_cov_hist) or force == True:
        filter_high_cover_profile(raw_hist, cds_fa, cover_ratio, cnt_threshold, high_cov_hist)
    else:
        print("high-coverage profile exists, use cached", file=sys.stderr)
    # step 3: train blur vector from meta profiles
    if not os.path.exists(vblur_fname) or force == True:
        train_vblur_from_meta(high_cov_hist, cds_fa, odir)
    else:
        print("vblur file exists, use cached", file=sys.stderr)
    # step 4: deblur high-coverage profiles
    if not os.path.exists(eps_fname) or force == True:
        deblur_transcripts(high_cov_hist, cds_fa, vblur_fname, odir)
    else:
        print("deblur file exits, use cached", file=sys.stderr)
    # step 5: combine length-specific profiles
    if not os.path.exists(profile_fname) or force == True:
        ctrue_merge = construct_deblur_profiles(eps_fname, high_cov_hist)
        i = 0
        for tid in ctrue_merge:
            print(tid, ctrue_merge[tid][:10])
            i += 1
            if i > 10: break
            
