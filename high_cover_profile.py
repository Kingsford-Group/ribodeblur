#!/usr/bin/env python
import sys
import argparse
import math
from global_params import *
from footprint_hist_parser import get_cds_range, parse_rlen_hist, get_transcript_profiles, write_rlen_hist

def get_high_cover_rlen_profile(plist, tlen, cover_ratio, cnt_threshold):
    """ 
    get high-coverage read count profiles per read lengh
    high-coverage profile: more than cover_ratio portion codon loci are high covered
    high-coverage location: read count > cnt_threshold
    codon loci: frame with the largest number of high-cover loci
    """
    transcript = {}
    # in-frame codon position coverage > cover_ratio
    pos_covered = math.ceil(tlen/3*cover_ratio)
    for rlen, profile in plist.items():
        pos_cnt = [0]*3
        for loc, cnt in profile:
            if cnt > cnt_threshold:
                pos_cnt[loc%3] += 1
        if max(pos_cnt) > pos_covered:
            transcript[rlen] = profile
    return transcript

def filter_transcript_profiles(tprofile, cds_range, cnt_threshold, cover_ratio):
    """
    filter high covered interesting profiles with different read lengths
    """
    print("filter transcript profiles", file=sys.stderr)
    pcelebrity = {}
    i = 0
    for tid, plist in tprofile.items():
        start, end = cds_range[tid]
        tlen = end-start
        transcript = get_high_cover_rlen_profile(plist, tlen, cover_ratio, cnt_threshold)
        if len(transcript)>1 and 28 in transcript:
            pcelebrity[tid] = transcript.copy()
        if 28 not in transcript and len(transcript)>0:
            print("WARNING:", tid, "no 28-mers", transcript.keys(), file=sys.stderr)
        i += 1
        print("processed transcript {}.\t\r".format(i), end="", file=sys.stderr, flush=True)
    print("\ntotal celebrity genes: {}".format(len(pcelebrity)), file=sys.stderr)
    return pcelebrity

def filter_high_cover_profile(hist_fn, cds_fa, cover_ratio, cnt_threshold, ofname):
    """ pipeline for filtering high coverage profiles """
    cds_range = get_cds_range(cds_fa)
    tlist = parse_rlen_hist(hist_fn)
    tprofile = get_transcript_profiles(tlist, cds_range, utr5_offset, utr3_offset)
    pcelebrity = filter_transcript_profiles(tprofile, cds_range, cnt_threshold, cover_ratio)
    tid2rid = { t['tid']: rid for rid, t in tlist.items() }
    write_rlen_hist(pcelebrity, cds_range, tid2rid, ofname)

def make_arg_parser():
    """ command line arguments """
    parser = argparse.ArgumentParser(prog="high_cover_profile", add_help=True, description="Filter and keep raw ribosome profiles with a high coverage")
    parser.add_argument("-i", "--input", required=True, help="Input rlen_hist file with all transcripts")
    parser.add_argument("-f", "--cds_fa", required=True, help="Reference fasta of padded transcriptome")
    parser.add_argument("-r", "--cover_ratio", type=float, default=0.5, help="High coverage ratio. Include transcript if it has more than this portion of loci with high coverage (default = 0.5)")
    parser.add_argument("-c", "--cnt_threshold", type=int, default=0, help="High coverage read count threshold. A location is high coverage if its read count is greater than such threshold (default = 0)")
    parser.add_argument("-o", "--output", required = True, help="Output rlen_hist file with high coverage transcripts")
    return parser

def main():
    """ command line wrapper """
    parser = make_arg_parser()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    filter_high_cover_profile(args.input, args.cds_fa, args.cover_ratio, args.cnt_threshold, args.output)

if __name__ == "__main__": main()
# python high_cover_profile.py -i test.hist -f ~/data/ribodeblur/refs/Saccharomyces_cerevisiae.R64-1-1.transcript_100.fa -o high_coverage.hist
