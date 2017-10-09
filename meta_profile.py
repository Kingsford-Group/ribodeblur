#!/usr/bin/env python
import numpy as np
import sys
from footprint_hist_parser import *
from global_params import *

#=============================
# utils
#=============================
def unzip_dic_tuples(dic):
    keys, vals = zip(*dic.items())
    order = np.argsort(keys)
    return np.array(keys)[order], np.array(vals)[order]

def get_max_len(cds_range):
    max_len = 0
    for tid, (start, end) in cds_range.items():
        if end > max_len:
            max_len = end
    return max_len
    
def is_sparse(vec, threshold):
    order = np.argsort(vec)
    # largest peak not included when determining sparsity
    if np.sum(vec[order][:-1])>= threshold:
        return False
    else:
        return True

def get_ylim(vec, pcnt=3):
    """ 
    pcnt: which peak to use as the ylim
    """
    order = np.argsort(vec)
    idx = order[-pcnt]
    return vec[idx]

def get_pos_hist(pos_hist, start, end):
    y_cnt = np.zeros(end-start+1)
    for pos, cnt in pos_hist.items():
        if pos < start or pos > end: continue
        idx = pos - start
        y_cnt[idx] = cnt
    return y_cnt

def get_cobs(pos_hist, rlstart, rlend, istart, iend):
    pobs = {}
    for rlen in range(rlstart, rlend+1):
        if rlen not in pos_hist: continue
        pobs[rlen] = get_pos_hist(pos_hist[rlen], istart, iend)
    return pobs

def get_frames(vec):
    fsum = np.array([ np.sum(vec[i::3]) for i in range(3) ])
    tsum = float(sum(fsum))
    f0 = np.argmax(fsum)
    fdic = { f0:0, (f0+1)%3:1, (f0+2)%3:2 }
    print("frames: {0:.2%} {1:.2%} {2:.2%}".format(fsum[0]/tsum, fsum[1]/tsum, fsum[2]/tsum))
    return fdic

def get_frame_str_vec(vec):
    fsum = np.array([ np.sum(vec[i::3]) for i in range(3) ])
    tsum = float(sum(fsum))
    f0 = np.argmax(fsum)
    fdic = { f0:0, (f0+1)%3:1, (f0+2)%3:2 }
    sframe =  "{0:.2%} {1:.2%} {2:.2%}".format(fsum[0]/tsum, fsum[1]/tsum, fsum[2]/tsum)
    return sframe  

def get_frame_str(pobs):
    return { rlen: get_frame_str_vec(vec) for rlen, vec in pobs.items() }


def get_frame_from_tlist(tlist, cds_range):
    sframe = {}
    for rid in tlist:
        for rlen in tlist[rid]['prof']:
            tid = tlist[rid]['tid']
            start, end = cds_range[tid]
            for pos, cnt in tlist[rid]['prof'][rlen]:
                sframe.setdefault(rlen, np.zeros(3))
                sframe[rlen][(pos-start)%3] += cnt
    return sframe
                
def get_frame_str_from_tlist(tlist, cds_range):
    sframe = get_frame_from_tlist(tlist, cds_range)
    return { rlen: "{0:.2%} {1:.2%} {2:.2%}".format(f[0]/sum(f), f[1]/sum(f), f[2]/sum(f)) for rlen,f in sframe.items() }

#=============================
# transcript pre-selection
#=============================
def filter_transcript_by_length(cds_range, length_cutoff):
    tid_list = np.array(list(cds_range.keys()))
    tlen = np.array([ cds_range[k][1]-cds_range[k][0] for k in tid_list])
    return set(tid_list[tlen > length_cutoff])

#=============================
# generate meta profile
#=============================
def create_rlen_meta_profile(tlist, cds_range, tid_included, ibegin, iend):
    """ 
    create count accumulation for each base location with different read length
    exclude transcript if id not in tid_included (a set)
    ibegin: index bound to include before START (negative)
    iend: number of bases to include after START (positive)
    """
    print("\ncreate meta hist")
    meta_hist = {}
    for rid in tlist:
        tid = tlist[rid]['tid']
        if tid not in tid_included: continue
        start, end = cds_range[tid]
        for rlen in tlist[rid]['prof']:
            meta_hist.setdefault(rlen, {})
            for pos, cnt in tlist[rid]['prof'][rlen]:
                i = pos - start
                if i < ibegin or i > iend : continue
                meta_hist[rlen].setdefault(i,0)
                meta_hist[rlen][i] += cnt
        sys.stdout.write("processed transcript {0}.\t\r".format(rid))
        sys.stdout.flush()
    sys.stdout.write("\n")
    return meta_hist

#=============================
# read length distribution
#=============================
def get_length_count(tlist):
    rlen2cnt = {}
    for rid in tlist:
        for rlen in tlist[rid]['prof']:
            tid = tlist[rid]['tid']
            for pos, cnt in tlist[rid]['prof'][rlen]:
                rlen2cnt.setdefault(rlen, 0)
                rlen2cnt[rlen] += cnt
    return rlen2cnt

def get_length_distribution(rlen2cnt, rlen_min, rlen_max):
    rlen2cnt_filter = { rlen:cnt for rlen,cnt in rlen2cnt.items() \
                        if rlen >= rlen_min and rlen <= rlen_max }
    tot_cnt = float(sum(rlen2cnt_filter.values()))
    return { rlen: cnt/tot_cnt for rlen,cnt in rlen2cnt_filter.items() }
