#!/usr/bin/env python
"""
build CDS with UTR padding from Ensembl genome fasta and gff annotation
prerequest python module: Biopython
"""
import argparse
import sys
import collections
import operator

import BCBio.GFF
import Bio.SeqIO
import Bio.SeqRecord

Exon = collections.namedtuple("Exon", ["start", "end", "strand", "phase", "phase_end", "rank"])

def extract_exon(seq_feature):
    """ seq_feature.type == exon, build exon info """
    start = int(seq_feature.location.start)
    end = int(seq_feature.location.end)
    strand = seq_feature.strand
    phase = int(seq_feature.qualifiers["ensembl_phase"][0])
    phase_end = int(seq_feature.qualifiers["ensembl_end_phase"][0])
    rank = int(seq_feature.qualifiers["rank"][0])
    return Exon(start, end, strand, phase, phase_end, rank)

def extract_parent(seq_feature):
    """ get parent (of exon feature) """
    tids =  seq_feature.qualifiers["Parent"]
    if len(tids) > 1:
        print("WARNING: more than one parent for exon {}! {}".format(seq_feature.id, tids))
    return tids[0]

def get_exon_list(seq_feature):
    """
    extract all exon children features of a seq feature
    """
    transcript = collections.defaultdict(list)
    # process queue for features
    feature_queue = collections.deque([seq_feature])
    while feature_queue:
        f = feature_queue.popleft()
        # feature is exon, append to transcript result
        if f.type == "exon":
            tid = extract_parent(f)
            exon = extract_exon(f)
            transcript[tid].append(exon)
        # feature is not exon, append sub_features to process queue
        else:
            feature_queue.extend(f.sub_features)
    return transcript

def get_annotations_from_gff(gff_fname):
    """ parse gff """
    gff_file = open(gff_fname)
    annotations = Bio.SeqIO.to_dict(BCBio.GFF.parse(gff_file))
    gff_file.close()
    return annotations

def extract_exons_from_annotations(annotations):
    """ 
    convert raw gff output into exon_list 
    output: { chrm: {transcript: [ Exon ] }}
    """
    transcriptome = collections.defaultdict(dict)
    for chrm in annotations:
        for feature in annotations[chrm].features:
            transcript = get_exon_list(feature)
            if transcript: 
                transcriptome[chrm].update(transcript)
    # sort exons by rank
    for chrm in transcriptome:
        for transcript in transcriptome[chrm]:
            exon_sorted = sorted(transcriptome[chrm][transcript],
                                 key = operator.attrgetter("rank"))
            transcriptome[chrm][transcript] = exon_sorted
    return transcriptome

def get_exons_from_gff(gff_fname):
    """ parse gff and output { chrm: {transcript: [ Exon ] }} """
    annotations = get_annotations_from_gff(gff_fname)
    transcriptome = extract_exons_from_annotations(annotations)
    return transcriptome

def append_padding(chrm_len, exon_list, padding):
    """ include UTR regions as 'Exons' """
    # adjust first exon
    e = exon_list[0]
    start = e.start - padding
    if start < 0: start = 0
    utr5 = Exon(start, e.start, e.strand, e.phase, e.phase_end, 0)
    # adjust last exon
    e = exon_list[-1]
    end = e.end + padding
    if end > chrm_len: end = chrm_len
    utr3 = Exon(e.end, end, e.strand, e.phase, e.phase_end, -1)
    return [utr5] + exon_list + [utr3]

def prepare_region_list(chrm_len, exon_list, padding):
    """
    reverse exon order if negative strand
    add UTR retions to exon_list
    """
    strand = exon_list[0].strand
    # reverse exon order for negative strand
    if strand == -1:
        exon_list = exon_list[::-1]
    # include UTR padding
    region_list = append_padding(chrm_len, exon_list, padding)
    return region_list

def get_cds(region_list):
    """ get CDS start and end from region_list """
    start = region_list[0].end - region_list[0].start
    end = region_list[-1].start - region_list[0].start
    return (start, end)

def get_seq(seq_rec, start, end):
    """ return regions of seq from seq_rec """
    return seq_rec.seq[start:end]

def connect_exons(seq_rec, exon_list):
    """
    return connected sequences from region list
    """
    seq = None
    for exon in exon_list:
        exon_seq = get_seq(seq_rec, exon.start, exon.end)
        if not seq: seq = exon_seq
        else: seq += exon_seq
    return seq

def build_transcript_seq(seq_rec, exon_list):
    """ return connected seq with strand handling """
    seq = connect_exons(seq_rec, exon_list)
    if exon_list[0].strand  == -1:
        return seq.reverse_complement()
    else: return seq

def build_transcriptome_rec(genome_seq, annotation, padding):
    """ build SeqRecord for transcripts """
    seq_rec_list = []
    for chrm in annotation:
        for transcript in annotation[chrm]:
            tid = transcript.lstrip("transcript:")
            region_list = prepare_region_list(len(genome_seq[chrm]), 
                                              annotation[chrm][transcript], 
                                              padding)
            seq = build_transcript_seq(genome_seq[chrm], region_list)
            start, end = get_cds(region_list)
            description = "|CDS:{}-{}|length:{}".format(start, end, len(seq))
            seq_rec_list.append(Bio.SeqRecord.SeqRecord(seq, 
                                                        id=tid, 
                                                        description=description))
    return seq_rec_list

def write_seq_to_fasta(seq_list, out_fa):
    """ write SeqRecords to fasta """
    ofile = open(out_fa, 'w')
    Bio.SeqIO.write(seq_list, ofile, "fasta")    
    ofile.close()
        
def make_arg_parser():
    """
    command line arguments
    """
    parser = argparse.ArgumentParser(prog="build_reference.py", add_help=True, \
                                     description="build coding region sequence fasta with utr padding from genome fasta and gff annotation")
    parser.add_argument("-f", "--input_fa", required=True, help="genome FASTA file from Ensembl")
    parser.add_argument("-a", "--gff", required=True, help="genome gff3 annotations from Ensembl")
    parser.add_argument("-o", "--output_fa", required=True, help="output FASTA of transcriptome")
    parser.add_argument("-p", "--padding", type=int, default=100, help="UTR length to be included outside of exons (default=100)")
    return parser

def main():
    """
    extract transcript sequences from genome fasta and gff annotations
    """
    print("preparing parameters...")
    parser = make_arg_parser()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    print("indexing genome fasta file...")
    genomes = Bio.SeqIO.index(args.input_fa, "fasta")
    print("parsing annotation gff...")
    exons = get_exons_from_gff(args.gff)
    print("getting transcriptome records...")
    seq_recs = build_transcriptome_rec(genomes, exons, args.padding)
    print("writing records to fasta...")
    write_seq_to_fasta(seq_recs, args.output_fa)    

if __name__ == "__main__": main()
