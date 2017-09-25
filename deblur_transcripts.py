#!/usr/bin/env python
import sys
import argparse
import collections
from multiprocessing import Pool
import numpy as np
from meta_profile import *
from deblur_utils import *
from deblur_result_io import *
from global_params import *
from map_to_reference import get_file_core

def batch_Asite_recovery(tprofile, cds_range, utr5_offset, utr3_offset, rlen_min, rlen_max, blur_vec, klist, converge_cutoff):
    ptrue = {}
    eps = {}
    i = 0
    for tid in tprofile:
        i += 1
        sys.stdout.write("deblur {0}th transcript: {1}\t\t\r".format(i, tid))
        sys.stdout.flush()        
        start, end = cds_range[tid]
        cobs = build_cobs_for_deblur(tprofile[tid], utr5_offset, (end-start)+utr3_offset, rlen_min, rlen_max)
        if len(cobs) == 0: continue
        ptrue_tid, eps_tid = recover_true_profile(cobs, klist, blur_vec, 0, 100, converge_cutoff)
        if np.all(ptrue_tid==0): continue
        ptrue[tid] = ptrue_tid
        eps[tid] = eps_tid
    print("\ntotal deblurred transcripts: {0}".format(len(ptrue)))
    return ptrue, eps

class single_transcript_asite(object):
    def __init__(self, blur_vec, klist, converge_cutoff):
        self.b = blur_vec
        self.k = klist
        self.c = converge_cutoff
    def __call__(self, params):
        tid = params[0]
        cobs = params[1]
        ptrue, eps = recover_true_profile(cobs, self.k, self.b, 0, 100, self.c)
        return tid, ptrue, eps

def batch_Asite_recovery_parallel(tprofile, cds_range, utr5_offset, utr3_offset, rlen_min, rlen_max, blur_vec, klist, converge_cutoff, nproc):
    cobs_all = np.array([ build_cobs_for_deblur(prof, utr5_offset, (cds_range[tid][1]-cds_range[tid][0])+utr3_offset, rlen_min, rlen_max) for tid,prof in tprofile.items() ])
    cobs_len_all = np.array(list(map(len, cobs_all)))
    tid_all = np.array(list(tprofile.keys()))
    cobs_in = cobs_all[cobs_len_all!=0]
    tid_in = tid_all[cobs_len_all!=0]
    pool = Pool(processes=nproc)
    results = [ r for r in pool.imap_unordered(single_transcript_asite(blur_vec, klist, converge_cutoff), zip(tid_in, cobs_in), 10) ]
    pool.close()
    pool.join()
    tid_list, ptrue_list, eps_list = zip(*results)
    tid_list = np.array(tid_list)
    ptrue_list = np.array(ptrue_list)
    eps_list = np.array(eps_list)
    valid = np.array(list(map(lambda x: x is not None, ptrue_list)))
    ptrue = dict(zip(tid_list[valid], ptrue_list[valid]))
    eps = dict(zip(tid_list[valid], eps_list[valid]))
    print("\ntotal deblurred transcripts: {0}".format(len(ptrue)))
    return ptrue, eps

def deblur_transcripts(hist_fn, cds_fa, vblur_fn, odir):
    ensure_dir(odir)
    print("get pre-computed blur vector")
    b = read_vblur(vblur_fn)
    vrlen_min, vrlen_max = get_rlen_range_from_vblur(b)
    cds_range = get_cds_range(cds_fa)
    tlist = parse_rlen_hist(hist_fn)
    # build profile for each transcript per read length
    tprofile = get_transcript_profiles(tlist, cds_range, utr5_offset, utr3_offset)
    print("batch A-site recovery")
    ptrue, eps = batch_Asite_recovery_parallel(tprofile, cds_range, utr5_offset, utr3_offset, vrlen_min, vrlen_max, b, klist, converge_cutoff, nproc)
    #ptrue, eps = batch_Asite_recovery(tprofile, cds_range, utr5_offset, utr3_offset, vrlen_min, vrlen_max, b, klist, converge_cutoff)
    ofname = odir+get_file_core(hist_fn)+".eps"
    write_essentials(ptrue, eps, ofname)

def make_arg_parser():
    """ command line arguments """
    parser = argparse.ArgumentParser(prog="deblur_transcripts.py", add_help=True, description="deblur length-specific profiles with high coverage")
    parser.add_argument("-i", "--input", required=True, help="Input rlen_hist file with high-coverage transcripts")
    parser.add_argument("-f", "--cds_fa", required=True, help="Reference fasta of padded transcriptome")
    parser.add_argument("-v", "--vblur", required=True, help="Blur vector parameters for deblur process")
    parser.add_argument("-o", "--output_dir", required = True, help="Output directoryto store blur vector")
    return parser

def main():
    """ parse bam file and group reads by length and transcript loci """
    parser = make_arg_parser()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    deblur_transcripts(args.input, args.cds_fa, args.vblur, args.output_dir)

if __name__ == "__main__": main()
