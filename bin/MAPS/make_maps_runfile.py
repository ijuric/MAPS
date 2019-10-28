##
## EXAMPLE RUN:
## python make_maps_runfile.py MY_113.MY_115_test /home/jurici/work/PLACseq/MAPS_pipe/results/mESC/ /home/jurici/work/PLACseq/analysis/binning_for_regression/new_set/final.replicated.narrowPeak /home/jurici/work/PLACseq/data/genomic_features/F_GC_M_MboI_5Kb_el.mm10.txt /home/jurici/work/PLACseq/analysis/binning_for_regression/new_set/ /home/jurici/work/PLACseq/analysis/binning_for_regression/new_set/ 5000 19 ./ NA 
##

import argparse
import sys

def init(parser):
    fname = str(parser.OUT_FILE_PATH) + 'maps_'+parser.DATASET_NAME+'.maps'
    outf = open(fname,'w')
    outstring = 'DATASET_NAME='+str(parser.DATASET_NAME)+'\n'+\
    'OUT_DIR='+str(parser.OUT_DIR)+'\n'+\
    'MACS2_PATH='+str(parser.MACS2_PATH)+'\n'\
    'GF_PATH='+str(parser.GF_PATH)+'\n'\
    'LONG_PATH='+str(parser.LONG_PATH)+'\n'\
    'LONG_FORMAT=[DATASET_NAME].[CHROMOSOME].long.intra.bedpe'+'\n'\
    'SHORT_PATH='+str(parser.SHORT_PATH)+'\n'\
    'SHORT_FORMAT=[DATASET_NAME].[CHROMOSOME].shrt.vip.bed'+'\n'\
    'BIN_SIZE='+str(parser.BIN_SIZE)+'\n'\
    'BINNING_RANGE='+str(parser.BINNING_RANGE)+'\n'\
    'N_CHROMS='+str(parser.N_CHROMS)+'\n'\
    'SEX_CHROMS='+str(parser.SEX_CHROMS)+'\n'
    outf.write(outstring)
    outf.close()

def main():
    parser = argparse.ArgumentParser()
    parser.prog = 'PROG'
    parser.description = "make_maps_runfile.py"
    parser.add_argument('DATASET_NAME', help = '')
    parser.add_argument('OUT_DIR', help = '')
    parser.add_argument('MACS2_PATH', help = '')
    parser.add_argument('GF_PATH', help = '')
    parser.add_argument('LONG_PATH', help = '')
    parser.add_argument('SHORT_PATH', help = '')
    parser.add_argument('BIN_SIZE', help = '')
    parser.add_argument('N_CHROMS', help = '')
    parser.add_argument('OUT_FILE_PATH', help = 'directory where .maps file will be stored. Must end with slash ')
    parser.add_argument('SEX_CHROMS', help='')
    parser.add_argument('--BINNING_RANGE', default=1000000, help = '')
    p = parser.parse_args(sys.argv[1:])
    init(p)

if __name__ == "__main__":
    main()
