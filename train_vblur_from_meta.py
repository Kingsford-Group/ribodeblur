#!/usr/bin/env python
import sys
import argparse
import numpy as np
from meta_profile import filter_transcript_by_length, create_rlen_meta_profile, get_cobs
from deblur_utils import train_vblur_from_meta_profiles
from deblur_result_io import write_vblur
from footprint_hist_parser import get_cds_range, parse_rlen_hist
from global_params import *

def get_max_frame_percent(vec, percentile=95):
    """ portion of reads in the most abundant frame """
    fsum = np.zeros(3)
    for i in range(3):
        vsub = vec[i::3]
        threshold = np.percentile(vsub, percentile)
        fsum[i] = np.sum(vsub[vsub<=threshold])
    if max(fsum)==0:
        return 0
    else:
        return max(fsum)/float(sum(fsum))

def get_min_mean_frame_cnt(vec, percentile=95):
    """ smallest average read count in 3 frames """
    fmean = np.zeros(3)
    for i in range(3):
        vsub = vec[i::3]
        threshold = np.percentile(vsub, percentile)
        fmean[i] = np.mean(vsub[vsub<=threshold])
    return min(fmean)

def get_vblur_rlen_range(mobs):
    """ 
    filter a range of read length to include
    profile of a given length need to be skew enough: most abundant frame > 40% within the 3 frames
    profile needs to have high coverage: least abundant frame > 10 read counts
    """
    vrlen_min = rlen_min
    vrlen_max = rlen_max
    for rlen in range(rlen_min, rlen_max+1):
        if rlen not in mobs: continue
        vrlen_min = rlen
        min_frame_cnt = get_min_mean_frame_cnt(mobs[rlen],100)
        max_frame_portion = get_max_frame_percent(mobs[rlen], 100)
        if min_frame_cnt >= lowest_frame_cnt and max_frame_portion >= lowest_frame_percent:
            break
    for rlen in range(rlen_max, rlen_min, -1):
        if rlen not in mobs: continue
        vrlen_max = rlen
        min_frame_cnt = get_min_mean_frame_cnt(mobs[rlen],100)
        max_frame_portion = get_max_frame_percent(mobs[rlen], 100)
        # print "{0} {1:.2%} {2:.0f}".format(rlen, max_frame_portion, min_frame_cnt)
        if min_frame_cnt >= lowest_frame_cnt and max_frame_portion >= lowest_frame_percent:
            break
    return vrlen_min, vrlen_max

def meta_pipeline(tlist, cds_range, istart, istop, rlen_min, rlen_max, converge_cutoff):
    """ train blur vector on meta profiles with different lengths """
    # train vblur on meta profile
    tid_select = filter_transcript_by_length(cds_range, istop)
    pos_hist = create_rlen_meta_profile(tlist, cds_range, tid_select, istart, istop)
    mobs = get_cobs(pos_hist, rlen_min, rlen_max, 0, istop)    
    vrlen_min, vrlen_max = get_vblur_rlen_range(mobs)
    mobs_hc = { rlen:mobs[rlen] for rlen in range(vrlen_min, vrlen_max+1) } 
    estep = True
    b, ptrue, eps = train_vblur_from_meta_profiles(mobs_hc, klist, low, percentile, converge_cutoff, estep)
    return b, ptrue, eps

def train_vblur_from_meta(hist_fn, cds_fa, vblur_fname):
    """ procedure pipeline """
    cds_range = get_cds_range(cds_fa)
    tlist = parse_rlen_hist(hist_fn)
    b, ptrue, eps = meta_pipeline(tlist, cds_range, utr5_offset, imax, rlen_min, rlen_max, converge_cutoff)
    write_vblur(b, vblur_fname)

def make_arg_parser():
    """ command line arguments """
    parser = argparse.ArgumentParser(prog="train_vblur_from_meta", add_help=True, description="Train blur vector from length-specific meta profiles")
    parser.add_argument("-i", "--input", required=True, help="Input rlen_hist file with high-coverage transcripts")
    parser.add_argument("-f", "--cds_fa", required=True, help="Reference fasta of padded transcriptome")
    parser.add_argument("-o", "--output", required = True, help="Output blur vector to file")
    return parser

def main():
    """ command line wrapper """
    parser = make_arg_parser()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    train_vblur_from_meta(args.input, args.cds_fa, args.output)

if __name__ == "__main__": main()
