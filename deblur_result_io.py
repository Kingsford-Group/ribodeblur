#!/usr/bin/env python
import os
import sys
import collections
import numpy as np

def ensure_dir(f):
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)

# def get_file_core(fname):
#     istart = fname.rfind("/")
#     iend = fname.rfind(".hist")
#     return fname[istart:iend]

def write_vblur(b, ofname):
    tf = open(ofname, 'w')
    text = [ "{0}: {1}\n".format(rlen, "\t".join(map(str,vblur)))
             for rlen, vblur in b.items() ]
    tf.writelines(text)
    tf.close()

def read_vblur(fname):
    tf = open(fname)
    b = collections.OrderedDict()
    for line in tf:
        read_len, pc_str = line.rstrip("\n").split(": ")
        read_len = int(read_len)
        prof = np.array(list(map(float, pc_str.split("\t"))))
        b[read_len] = prof
    tf.close()
    return b
    
def write_essentials(ptrue, eps, ofname):
    if len(ptrue)==0:
        return
    tf = open(ofname, 'w')
    i = 0
    for tid, prof in ptrue.items():
        text = [ "tid: {0}\n".format(tid), 
                 "0: {0}\n".format("\t".join(map(str, prof))) ] + \
        [ "{0}: {1}\n".format(rlen, "\t".join(map(str, eps_rlen)))
          for rlen, eps_rlen in eps[tid].items() ]
        tf.writelines(text)
        i += 1
        sys.stdout.write("processed transcripts {0}.\t\r".format(i))
        sys.stdout.flush()
    sys.stdout.write("\n")
    tf.close()

def read_essentials(fname):
    tf = open(fname)
    ptrue = {}
    eps = {}
    transcript = {}
    line = tf.readline()
    i = 0
    while line:
        if line.startswith("tid: "):
            if transcript:
                eps[tid] = transcript.copy()
                i += 1
                sys.stdout.write("processed transcripts {0}.\t\r".format(i))
                sys.stdout.flush()
                transcript.clear()
            tid = line.lstrip("tid: ").rstrip("\n")
        else:
            read_len, pc_str = line.rstrip("\n").split(": ")
            read_len = int(read_len)
            prof = np.array(list(map(float, pc_str.split("\t"))))
            if read_len == 0:
                ptrue[tid] = prof
            else:
                transcript[read_len] = prof
        line = tf.readline()
    if transcript:
        eps[tid] = transcript.copy()
    sys.stdout.write("\n")
    tf.close()
    return ptrue, eps

def write_cds_profile(base_prof, ofname):
    if len(base_prof) == 0:
        return
    tf = open(ofname, 'w')
    i = 0
    for tid, prof in base_prof.items():
        text = [ "tid: {0}\n".format(tid), 
                 "{0}\n".format("\t".join(map(str, prof))) ]
        tf.writelines(text)
        i += 1
        sys.stdout.write("processed transcripts {0}.\t\r".format(i))
        sys.stdout.flush()
    sys.stdout.write("\n")
    tf.close()

def read_cds_profile(fname):
    tf = open(fname)
    base_prof = {}
    line = tf.readline()
    i = 0
    while line:
        if line.startswith("tid: "):
            i += 1
            tid = line.lstrip("tid: ").rstrip("\n")
            line = tf.readline()
            prof = map(float, line.rstrip("\n").split("\t"))
            base_prof[tid] = np.array(prof)
            sys.stdout.write("processed transcripts {0}.\t\r".format(i))
            sys.stdout.flush()
        line = tf.readline()
    sys.stdout.write("\n")
    tf.close()
    return base_prof

#=============================
# format conversion
#=============================
def build_cobs_per_rlen_with_shifts(pos_list, start, end, offset):
    """ shift left with offset"""
    profile = np.zeros(end-start+1)
    for pos, cnt in pos_list:
        idx_adj = pos-start-offset
        if idx_adj < 0 or idx_adj >= len(profile): continue
        profile[idx_adj] = cnt
    return profile

def build_cobs_per_transcript_with_shifts(clist, start, end, rlen_min, rlen_max, klist):
    cobst = {}
    for rlen, pos_list in clist.items():
        if rlen < rlen_min or rlen > rlen_max: continue
        cobst[rlen] = build_cobs_per_rlen_with_shifts(pos_list, start, end, klist[rlen])
    return cobst

def build_cobs_with_shifts(tprofile, cds_range, utr5_offset, utr3_offset, rlen_min, rlen_max, klist):
    cobs = {}
    i = 0
    for tid, prof in tprofile.items():
        start, end = cds_range[tid]
        cobs[tid] = build_cobs_per_transcript_with_shifts(prof, utr5_offset, (end-start)+utr3_offset, rlen_min, rlen_max, klist)
        i += 1
        sys.stdout.write("processed transcript {0}.\t\r".format(i))
        sys.stdout.flush()
    sys.stdout.write("\n")
    return cobs

def build_profile_from_list(pos_list, start, end):
    profile = np.zeros(end-start+1)
    for pos, cnt in pos_list:
        profile[pos-start] = cnt
    return profile

def build_cobs_for_deblur(clist, start, end, rlen_min, rlen_max):
    cobs = {}
    for rlen, pos_list in clist.items():
        if rlen < rlen_min or rlen > rlen_max: continue
        cobs[rlen] = build_profile_from_list(pos_list, start, end)
    return cobs
