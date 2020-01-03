#!/usr/bin/python

#########################################
# Author: Yunjiang Qiu <serein927@gmail.com>
# File: merge_map.py
# Create Date: 2015-10-09 12:27:31
#########################################

import sys
import argparse

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-c", "--cut", dest="cut",required=True,help="RE cut file")
    parser.add_argument("-m", "--map", dest="map",required=True,help="mapability file")
    parser.add_argument("-o", "--output", dest="outfile",required=True,help="output file")
    args = parser.parse_args()

    mapability = {}
    with open(args.map,'r') as f:
        for line in f:
            tmp = line.rstrip().split('\t')
            mapability[tmp[0]] = tmp[4]

    with open(args.outfile,'w') as fo:
        with open(args.cut,'r') as f:
            for line in f:
                tmp = line.rstrip().split('\t')
                key = tmp[2] + "_" + tmp[0]
                try:
                    fo.write("%s\n" % "\t".join([str(i) for i in tmp + [mapability[key]]]))
                except KeyError:
                    sys.exit("Frag %s cannot be found in mapability file." % key)

if __name__ == "__main__":
    sys.exit(main())

