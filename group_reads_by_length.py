#!/usr/bin/env python
"""
group reads by read length
stored as a list of read start position and read count per transcript
"""

import sys
import argparse
import collections
import pysam

# cigar string operator map
# from: http://pysam.readthedocs.org/en/latest/api.html#pysam.AlignedSegment
# M BAM_CMATCH  0
# I BAM_CINS    1
# D BAM_CDEL    2
# N BAM_CREF_SKIP   3
# S BAM_CSOFT_CLIP  4
# H BAM_CHARD_CLIP  5
# P BAM_CPAD    6
# = BAM_CEQUAL  7
# X BAM_CDIFF   8
cop = {op:i for i,op in enumerate("MIDNSHP=X")}

def get_read_len(read_rec):
    """ get read length from parsing cigar string """
    read_len = 0
    for op,l in read_rec.cigar:
        # splice: indicate novel splices, skip
        if op == 3: return 0
        # only match/mismatch or deletion will change read_len
        if op == 0 or op == 2:
            read_len += l
    return read_len

def is_multimapped(read_rec):
    """ return True if multi-mapped ("NH" > 1) """
    if read_rec.has_tag("NH"):
        if read_rec.get_tag("NH") > 1: 
            return True
    return False

class ReadLenHist():
    """ read counts grouped by transcript, read length, and positions """
    def __init__(self, bam_obj):
        """ get transcript name and initialize counter """
        self.ref_names = bam_obj.references
        # reference ID : { rlen : Counter }
        self.hist = [ collections.defaultdict(collections.Counter) 
                      for i in range(bam_obj.nreferences) ]
        self.nref = bam_obj.nreferences
        self.mapped = 0
        self.unmapped = 0

    def update_counter(self, read):
        """ update read count given an alignment from pysam """
        # ignore reverse complement and multimapping
        if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_reverse or is_multimapped(read): 
            self.unmapped += 1
        else:
            rlen = get_read_len(read)
            if rlen == 0: 
                self.unmapped += 1
            else:
                # read start position on transcript
                # 0-based left-most coordinate
                self.mapped += 1
                self.hist[read.reference_id][rlen].update([read.reference_start])

    def summarize_reads(self):
        """ summarize mapping statistics of a read set """
        self.tot = self.mapped + self.unmapped
        print("\ntotal: {:,} uniquely mapped {:,}({:.2%}) rest: {:,} ({:.2%})".format(self.tot, self.mapped, self.mapped/self.tot, self.unmapped, self.unmapped/self.tot), file=sys.stderr)

def bam_to_rlen_group(bam_fname):
    """ parse bam file and group reads by length """
    bam = pysam.AlignmentFile(bam_fname, 'rb')
    rlen_hist = ReadLenHist(bam)
    i = 0
    for read in bam.fetch(until_eof=True):
        i += 1
        if i % 100000 == 0:
            print("\r\rprocessed {0:,} reads...".format(i), 
                  end=" ", file=sys.stderr, flush=True)
        rlen_hist.update_counter(read)
    rlen_hist.summarize_reads()
    return rlen_hist

def write_rlen_hist(rlen_hist, ofname):
    """ write read count grouped by length to file """
    ofile = open(ofname, 'w')
    text = []
    for i in range(rlen_hist.nref):
        if rlen_hist.hist[i]:
            text.append("refID: {}".format(i))
            text.append("tid: {}".format(rlen_hist.ref_names[i]))
            for rlen in sorted(rlen_hist.hist[i].keys()):
                pos_cnt_list = ["{} {}".format(k,rlen_hist.hist[i][rlen][k])
                                for k in sorted(rlen_hist.hist[i][rlen].keys())]
                text.append("{}: {}".format(rlen, ", ".join(pos_cnt_list)))
    ofile.write("\n".join(text))
    ofile.close()

def group_reads_by_length(ifname, ofname):
    """ pipeline for grouping reads by length """
    rlen_hist = bam_to_rlen_group(args.bam)
    write_rlen_hist(rlen_hist, args.ofname)
    
def make_arg_parser():
    """ command line arguments """
    parser = argparse.ArgumentParser(prog="group_reads_by_length.py", add_help=True, description="group read alignments by read length and transcript positions")
    parser.add_argument("-b", "--bam", required=True, help="read alignment bam file")
    parser.add_argument("-o", "--ofname", required = True, help="output file name to store read counts grouped by read length and transcript positions")
    return parser

def main():
    """ parse bam file and group reads by length and transcript loci """
    parser = make_arg_parser()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    group_reads_by_len_pipeline(args.bam, args.ofname)

if __name__ == "__main__": main()
# python group_reads_by_length.py -b /home/hw1/data/ribodeblur/star_align/SRR1177157_transcript_Aligned.out.bam -o raw.hist
