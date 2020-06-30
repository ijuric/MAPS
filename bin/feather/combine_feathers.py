import argparse
import os
import glob
import re
import subprocess

def main():
    parser = create_args()
    args = parser.parse_args()
    chroms = get_union_chroms(args.datasets)
    print('combining files from', len(chroms), 'chromosomes over', len(args.datasets), 'datasets...')
    combine_longs(args.datasets, chroms, args.outdir, args.prefix)
    combine_shorts(args.datasets, chroms, args.outdir, args.prefix)
    print('computing QC metrics...')
    compute_qc(args.datasets, args.outdir, args.prefix, args.chipseq)

def compute_qc(datasets, outdir, prefix, chipseq):
    int_keys = ['fastq_all', 'mapped', 'not_duplicated',             'intra_all', 'short_all', 'short_final', 'long_intra', 'inter_all']
    ratio_keys = ['trans_ratio', 'long_cis_ratio', 'FRiP']
    qc = {}
    for i, dataset in enumerate(datasets):
        qc[i] = {}
        qc_file = glob.glob(f'{dataset}/*.feather.qc.tsv')[0]
        with open(qc_file) as infile:
            lines = infile.readlines()
        for line in lines:
            line = line.split()
            qc[i][line[0]] = [line[1], line[2]]
    combined_qc = {}
    for key in int_keys:
        combined_qc[key] = sum([int(float(qc[ds][key][0])) for ds in qc.keys()])
    combined_qc['trans_ratio'] = "{0:.5f}".format(float(combined_qc['inter_all']) / combined_qc['not_duplicated'])
    combined_qc['long_cis_ratio'] = "{0:.5f}".format(float(combined_qc['long_intra']) / combined_qc['intra_all'])
    combined_qc['FRiP'] = "{0:.5f}".format(compute_frip(f'{outdir}/{prefix}.shrt.vip.bed', chipseq, outdir, prefix)) if chipseq else -1
    with open(f'{outdir}/{prefix}.feather.qc.tsv', 'w') as ofile:
        example_qc = qc[list(qc.keys())[0]]
        for keyname in int_keys + ratio_keys:
            ofile.write('\t'.join([keyname, str(combined_qc[keyname]), example_qc[keyname][1]]) + '\n')

def compute_frip(short_filename, chip_filename, outdir, prefix):
    cut_short_filename = short_filename + ".cut"
    command = "cut -f1-3 " + short_filename + " > " + cut_short_filename
    proc = subprocess.Popen(command, stdout = subprocess.PIPE, shell= True)
    proc.communicate()
    cut_chip_filename = f'{outdir}/chip_peaks.cut'
    command = "cut -f1-3 " + chip_filename + " > " + cut_chip_filename
    proc = subprocess.Popen(command, stdout = subprocess.PIPE, shell= True)
    proc.communicate()
    command = "bedtools intersect -a " + cut_short_filename + " -b " + cut_chip_filename + " | wc -l"
    proc = subprocess.Popen(command, stdout = subprocess.PIPE, shell= True)
    on_chip = int(proc.stdout.read().decode("utf-8"))
    command = "wc -l " + cut_short_filename
    #print (command)
    proc = subprocess.Popen(command, shell = True, stdout = subprocess.PIPE)
    chip_count = int(proc.stdout.read().decode("utf-8").split()[0])
    #print(on_chip, chip_count)
    on_chip_ratio = float(on_chip) / chip_count 
    return on_chip_ratio

def extract_filenames(datasets, chrom, suffix):
    filenames = [glob.glob(f'{dataset}/*.{chrom}.{suffix}') for dataset in datasets]
    file_per_datasets = list(map(len, filenames))
    if max(file_per_datasets) > 1:
        raise Exception(f'inconsistent naming problem ({suffix})')
    filenames = [name[0] for name in filenames if len(name) == 1]
    return filenames
    
def cat_files(filenames, outfile):
    cmd = 'cat ' + ' '.join(filenames) + ' > ' + outfile
    proc = subprocess.Popen(cmd, shell = True)
    proc.communicate()

def combine_longs(datasets, chroms, outdir, prefix):
    output_prefix = f'{outdir}/{prefix}'
    for chrom in chroms:
        outfile = f'{output_prefix}.{chrom}.long.intra.bedpe'
        filenames = extract_filenames(datasets, chrom, 'long.intra.bedpe')
        cat_files(filenames, outfile)

def combine_shorts(datasets, chroms, outdir, prefix):
    output_prefix = f'{outdir}/{prefix}'
    shrt_files = []
    for chrom in chroms:
        outfile = f'{output_prefix}.{chrom}.shrt.vip.bed'
        filenames = extract_filenames(datasets, chrom, 'shrt.vip.bed')
        cat_files(filenames, outfile)
        shrt_files.append(outfile)
    cat_files(shrt_files, f'{output_prefix}.shrt.vip.bed')
    
def get_union_chroms(datasets):
    chroms = set()
    for dataset in datasets:
        #print('looking for ' + f'{dataset}/*.chr*.long.intra.bedpe')
        filenames = glob.glob(f'{dataset}/*.chr*.long.intra.bedpe')
        chroms = chroms.union(set(map(lambda x: re.search('chr.*?\.', x)[0][:-1], filenames)))
    return chroms
    
def create_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--prefix', help = 'dataset name')
    parser.add_argument('-o', '--outdir', help = 'output directory')
    parser.add_argument('-d', '--datasets', nargs = '+', help = 'space separated paths to directories to be merged')
    parser.add_argument('-a', '--chipseq', help = 'filepath for ChIP-seq bed file, used for QC', required = False, default = None)
    return parser

if __name__ == '__main__':
    main()
