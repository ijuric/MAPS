#!/usr/bin/python

#########################################
# Author: Yunjiang Qiu <serein927@gmail.com>
# File: feature_frag2bin.py
# Create Date: 2015-07-07 18:11:28
#########################################

import sys
import argparse

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-i", "--input", dest="infile",required=True,help="input file")
    parser.add_argument("-b", "--bin_size", dest="bin_size",required=True,help="bin_size")
    parser.add_argument("-g", "--genome_size", dest="g_size",required=True,help="genome_size")
    parser.add_argument("-o", "--output", dest="outfile",required=True,help="output file")
    args = parser.parse_args()

    bin_size = args.bin_size.replace('Kb','000')
    bin_size = bin_size.replace('Mb','000000')
    try:
        bin_size = int(bin_size)
    except ValueError:
        sys.exit("Unknown bin size %s, please double check." % args.bin_size)

    g_size = {}
    with open(args.g_size, 'r') as f:
        for line in f:
            key,value = line.rstrip().split('\t')
            g_size[key] = int(value)

    feature_frag = {}
    with open(args.infile,'r') as f:
        for line in f:
            feat = line.rstrip().split('\t')
            assert feat[2] in g_size,"%s is not in genome size files!" % feat[2]
            assert int(feat[3]) <= g_size[feat[2]],"%d is larger than size of %s!" % (int(feat[3]),feat[2])
            key = '\t'.join([feat[2],str((int(feat[3])//bin_size))])
            value = (int(feat[4]),float(feat[5]),float(feat[6])) #frag length,GC,mappability
            try:
                feature_frag[key].append(value)
            except KeyError:
                feature_frag[key] = [value]
    #print(list(feature_frag.keys())[1:5])
    with open(args.outfile,'w') as f:
        for chr_name in sorted(g_size.keys(),key=lambda i:i[2:]):
            for bin_num in range(int(g_size[chr_name]/bin_size+1)):
                key = '\t'.join([chr_name,str(bin_num)])
                try:
                    #print(key)
                    xtmp = feature_frag['\t'.join([chr_name,str(bin_num)])]
                    #print("found")
                    #print(len(xtmp))
                    xtmplen = list(set([i[0] for i in xtmp]))
                    gc = sum([i[1] for i in xtmp]) / len(xtmp)
                    frag_len = sum([i for i in xtmplen if i<1000]) + len([i for i in xtmplen if i>=1000]) * 1000
                    mappability = sum([i[2] for i in xtmp]) / len(xtmp)
                except KeyError:
                    #print("not found")
                    gc = 0
                    frag_len = 0
                    mappability = 0
                f.write('%s\t%d\t%d\t%d\t%.4f\t%.4f\n' % (chr_name,bin_num*bin_size,(bin_num+1)*bin_size,frag_len,gc,mappability))

if __name__ == "__main__":
    sys.exit(main())

