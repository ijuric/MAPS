#!/usr/bin/env python

"""
split bam file into .shrt.bam and .long.bam
Created by Rongxin Fang
Modified by Armen Abnousi
"""
import sys
import pysam 
import collections
from itertools import islice
import pybedtools
import time
import sys

def bam2bed(bam1, bed1):
	print(time.ctime() + "starting the splitting operation")
	long_inter_bam = pybedtools.BedTool(bam1)
	long_inter_bedpe = long_inter_bam.bam_to_bed(bedpe = True, stream = True)
	long_inter_bedpe.saveas(bed1)
	print(time.ctime() + " splitting completed")
bam1 = sys.argv[1]
bed1 = sys.argv[2]
print("convering " + bam1)
bam2bed(bam1, bed1)
