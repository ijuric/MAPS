#!/usr/bin/python

#########################################
# Author: Yunjiang Qiu <serein927@gmail.com>
# File: restriction_cut.py
# Create Date: 2015-02-27 16:46:15
#########################################

import sys
import argparse
from Bio import SeqIO
from Bio.SeqUtils import GC
import re

def find_site(fasta,seq,outfile,pos):
    f = open(outfile, 'w')
    size = len(seq)
    for seq_record in SeqIO.parse(fasta, "fasta"):
        sys.stderr.write("processing "+seq_record.id+"\n")
        site = [0 - pos] + [m.start() + 1 for m in re.finditer(seq.lower(), str(seq_record.seq).lower())] + [len(str(seq_record.seq)) + 1 + pos - size]

        for i in range(1,len(site)-1):
            count = (i -1) * 2 + 1
            frag_start = site[i-1] + pos
            frag_end = site[i] + size - pos - 1
            frag_len = frag_end - frag_start
            frag_gc = GC(str(seq_record.seq)[max(frag_end-200,0):frag_end]) / 100
            f.write("{num}\t{strand}\t{chr}\t{pos}\t{fraglen}\t{GC}\n".format(num=count,strand="-",chr=seq_record.id,pos=frag_end,fraglen=frag_len,GC=frag_gc))
            frag_start = site[i] + pos
            frag_end = site[i+1] + size - pos - 1
            frag_len = frag_end - frag_start
            frag_gc = GC(str(seq_record.seq)[frag_start:frag_start+200]) / 100
            f.write("{num}\t{strand}\t{chr}\t{pos}\t{fraglen}\t{GC}\n".format(num=count+1,strand="+",chr=seq_record.id,pos=frag_start,fraglen=frag_len,GC=frag_gc))

    f.close()

def main():
    parser = argparse.ArgumentParser(description='')
    parser = argparse.ArgumentParser(description='Given a fasta file and a restricted recognition site, generate genomic features')
    parser.add_argument("-f", "--fasta", dest="fasta",required=True,help="input fasta file")
    parser.add_argument("-s", "--seq", dest="seq",required=True,help="RE cut sequence")
    parser.add_argument("-o", "--out", dest="outfile",required=True,help="Output file")
    parser.add_argument("-p", "--pos", dest="pos",required=True,help="RE cut position")
    args = parser.parse_args()
    find_site(args.fasta,args.seq,args.outfile,int(args.pos))

if __name__ == "__main__":
    sys.exit(main())
