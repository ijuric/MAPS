#!/usr/bin/python

#########################################
# Author: Yunjiang Qiu <serein927@gmail.com>
# File: restriction_cut.py
# Create Date: 2015-02-27 16:46:15
#########################################

import sys
import argparse
from Bio import SeqIO
import re

def find_site(fasta,seq,outfile):
    f = open(outfile, 'w')
    for seq_record in SeqIO.parse(fasta, "fasta"):
        sys.stderr.write("processing "+seq_record.id+"\n")
        site = [m.start() + 1 for m in re.finditer(seq.lower(), str(seq_record.seq).lower())]
        f.write("%s " % seq_record.id.replace("chr",""))
        for i in site:
            f.write("%s " % i)

        f.write("%d\n" % len(str(seq_record.seq)))
    f.close()

def main():
    parser = argparse.ArgumentParser(description='')
    parser = argparse.ArgumentParser(description='Given a fasta file and a restricted recognition site, generate genomic features')
    parser.add_argument("-f", "--fasta", dest="fasta",required=True,help="input fasta file")
    parser.add_argument("-s", "--seq", dest="seq",required=True,help="RE cut sequence")
    parser.add_argument("-o", "--out", dest="outfile",required=True,help="Output file")
    args = parser.parse_args()
    find_site(args.fasta,args.seq,args.outfile)

if __name__ == "__main__":
    sys.exit(main())
