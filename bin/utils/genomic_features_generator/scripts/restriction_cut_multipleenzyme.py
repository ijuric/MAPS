#!/usr/bin/python

#########################################
# Author: Yunjiang Qiu <serein927@gmail.com>
# Modified by: Armen Abnousi
# File: restriction_cut.py
# Create Date: 2015-02-27 16:46:15
#########################################

import sys
import argparse
from Bio import SeqIO
from Bio.SeqUtils import GC
import re

def find_site(fasta,seq,outfile,pos, binsize):
    f = open(outfile, 'w')
    mnase = False
    if (len(seq.split(",")) > 1):
        seqs = seq.split(",")
        sizes = [len(seq) for seq in seqs]
        poses = pos.split(",")
        poses = [int(x) for x in poses]
    else:
        mnase = True if seq == "mnase" else False
        seqs = [seq]
        sizes = [len(seq)]
        poses = [int(pos)] if not mnase else [0]
    site = []
    site_pos = []
    site_size = []
    for seq_record in SeqIO.parse(fasta, "fasta"):
        #sys.stderr.write("processing "+ str(mnase)+"\n")
        #sys.stderr.write("processing "+seq_record.id+"\n")
        #sys.stderr.write("processing "+str(len(seq_record.seq))+" "+str(binsize)+"\n")
        if not mnase:
            sys.stderr.write("processing "+seq_record.id+"\n")
            #print(seq_record.seq[3000181: (3000185 + 30)])
            site = []
            site_pos = []
            site_size = []
            #so [0-pos] and [len(str(seq_record.seq)) + 1 + pos - size] are added so we can work with them in the for loop, instead of [0-pos] do [0-pos for the first seq found in the seq_record] and instead of pos and size in the second one use pos and size specific to the last re found in the seq_record
            for i in range(len(seqs)):
                seq = seqs[i]
                size = sizes[i]
                pos = poses[i]
                #print("searching for:" + seq.lower())
	        ###this is what needs to change:
                current_site = [m.start() + 1 for m in re.finditer(seq.lower(), str(seq_record.seq).lower())]
                site = site + current_site
                site_pos = site_pos + ([pos] * len(current_site))
                site_size = site_size + ([size] * len(current_site))
                ##for each site record what is its size and pos and then zip and sort by site, you need the site and pos in the loop below
            if (len(site) == 0):
                continue
            site, site_pos, site_size = (list(x) for x in zip(*sorted(zip(site, site_pos, site_size))))
            #print(site[:10])
            #break
            first_pos = site_pos[0]
            first_size = site_size[0]
            last_pos = site_pos[-1]
            last_size = site_size[-1]
            site = [0 - first_pos] + site + [len(str(seq_record.seq)) + 1 + last_pos - last_size]
            site_pos = [first_pos] + site_pos + [last_pos]
            site_size = [first_size] + site_size + [last_size]

            for i in range(1,len(site)-1):
                count = (i -1) * 2 + 1
                frag_start = site[i-1] + site_pos[i-1]
                frag_end = site[i] + site_size[i] - site_pos[i] - 1
                frag_len = frag_end - frag_start
                frag_gc = GC(str(seq_record.seq)[max(frag_end-200,0):frag_end]) / 100
                f.write("{num}\t{strand}\t{chr}\t{pos}\t{fraglen}\t{GC}\n".format(num=count,strand="-",chr=seq_record.id,pos=frag_end,fraglen=frag_len,GC=frag_gc))
                frag_start = site[i] + site_pos[i]
                frag_end = site[i+1] + site_size[i+1] - site_pos[i+1] - 1
                frag_len = frag_end - frag_start
                frag_gc = GC(str(seq_record.seq)[frag_start:frag_start+200]) / 100
                f.write("{num}\t{strand}\t{chr}\t{pos}\t{fraglen}\t{GC}\n".format(num=count+1,strand="+",chr=seq_record.id,pos=frag_start,fraglen=frag_len,GC=frag_gc))
        elif mnase:
            for i, frag_start in enumerate(range(1, len(str(seq_record.seq)) + 1, binsize)):
                count = i * 2 + 1
                frag_end = min(frag_start + binsize - 1, len(str(seq_record.seq)))
                frag_len = binsize - 1
                frag_gc = GC(str(seq_record.seq)[frag_start:frag_end]) / 100
                f.write("{num}\t{strand}\t{chr}\t{pos}\t{fraglen}\t{GC}\n".format(num=count,strand="-",chr=seq_record.id,pos=frag_end,fraglen=frag_len,GC=frag_gc))
                f.write("{num}\t{strand}\t{chr}\t{pos}\t{fraglen}\t{GC}\n".format(num=count+1,strand="+",chr=seq_record.id,pos=frag_start,fraglen=frag_len,GC=frag_gc))

    f.close()

def main():
    parser = argparse.ArgumentParser(description='')
    parser = argparse.ArgumentParser(description='Given a fasta file and a restricted recognition site, generate genomic features')
    parser.add_argument("-f", "--fasta", dest="fasta",required=True,help="input fasta file")
    parser.add_argument("-s", "--seq", dest="seq",required=True,help="RE cut sequence")
    parser.add_argument("-o", "--out", dest="outfile",required=True,help="Output file")
    parser.add_argument("-p", "--pos", dest="pos",required=True,help="RE cut position")
    parser.add_argument("-b", "--binsize", dest = "binsize", required = False, help = "bin size for MNase-based", default = "5Kb")
    args = parser.parse_args()
    bin_size = args.binsize.replace('Kb','000')
    bin_size = bin_size.replace('Mb','000000')
    try:
        bin_size = int(bin_size)
    except ValueError:
        sys.exit("Unknown bin size %s, please double check." % args.bin_size)
    find_site(args.fasta,args.seq,args.outfile,args.pos, bin_size)

if __name__ == "__main__":
    sys.exit(main())
