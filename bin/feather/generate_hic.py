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
from itertools import combinations
import time
import os
import subprocess
import shutil
import glob

def main():
	input_bam = sys.argv[1]
	output = sys.argv[2]
	setname = sys.argv[3]
	outdir = output[:max(0, output.rfind("/"))]
	print("\noutput: " + output + "\n")
	print("\noutdir: " + outdir + "\n")
	if (len(outdir) == 0):
		outdir = "."
	split_main(input_bam, outdir, setname)
	temp_files = glob.glob(outdir + "/hic_tempfiles_" + setname + "/*")
	command = ["cat"]
	command.extend(temp_files)
	print("catting per chromosome pair files into a single file.")
	with open(output, 'w') as outfile:
		proc = subprocess.Popen(command, stdout = outfile, shell= False)
		proc.communicate()

def split_main(input_bam, outdir, setname):
	print(time.ctime() + " starting the splitting operation")
	print(input_bam)
	samfile = pysam.AlignmentFile(input_bam, "rb")
	chr_list = extract_chr_list(samfile.header)
	if (len(chr_list) > 40):
		print("NOTE: Too many chromosomes. Generating hic only for ones without an underscore in their name")
		chr_list = [chr_name for chr_name in chr_list if (chr_name.find('_') == -1)]
		#print("\n".join(chr_list))
	if not is_sorted_queryname(samfile.header):
		sys.exit("Error: bam needs to be sorted by read name")
	hic_tsv_files = open_hic_files(chr_list, outdir, setname)
	read_num = 0
	prev = pysam.AlignedSegment()
	for read in samfile:
		if read_num % 2 == 1:
			if read.query_name != prev.query_name:
				exit("unmatched query" + read.query_name)
		if (prev.reference_name not in chr_list) or (read.reference_name not in chr_list):
			prev = read
			read_num += 1
			continue
		pos1 = prev.reference_end if prev.is_reverse else prev.reference_start
		pos2 = read.reference_end if read.is_reverse else read.reference_start
		if read_num % 2 == 1:
			'''
			try:
				chrom_num1 = int(float(read.reference_name[3:]))
				chrom_num2 = int(float(prev.reference_name[3:]))
				if chrom_num1 > chrom_num2:
					r1 = prev
					r2 = read
				elif chrom_num1 == chrom_num2 and read.reference_start > prev.reference_start:
					r1 = prev
					r2 = read
				else:
					r1 = read
					r2 = prev
			except ValueError:
				r1 = read
				r2 = prev
			'''
			if read.reference_start > prev.reference_start:
				r1 = prev
				r2 = read
			else:
				r2 = prev
				r1 = read
			chr1 = r1.reference_name
			chr2 = r2.reference_name
			if (r1.is_reverse and r2.is_reverse):
				strand1 = 16
				strand2 = 16
			elif (r1.is_reverse and r2.is_reverse == False):
				strand1 = 16
				strand2 = 0
			elif (r1.is_reverse == False and r2.is_reverse):
				strand1 = 0
				strand2 = 16
			else:
				strand1 = 0
				strand2 = 0
			hic_tsv_files[tuple(sorted([chr1, chr2]))].write("\t".join([r1.query_name, str(strand1), chr1, str(r1.reference_start), "0", str(strand2), chr2, str(r2.reference_start), "1", "60", "60\n"]))
		read_num += 1
		prev = read
	close_files(hic_tsv_files)
	print(time.ctime() + " per chromosome pair hic input generation completed")

def close_files(files_dict):
	for filename in files_dict.values():
		filename.close()

def open_hic_files(chr_list, outdir, setname):
	files_dir =  outdir + "/hic_tempfiles_" + setname
	print("\ncreating dir: " + files_dir + "\n")
	#print(len(chr_list))
	if os.path.exists(files_dir):
		shutil.rmtree(files_dir)
	os.mkdir(files_dir, 0o777)
	file_prefix = files_dir + "/hic."
	files_dict = {}
	for i, j in combinations(chr_list, 2):
		i, j = sorted([i, j])
		fout1 = open(file_prefix + i + "." + j, 'w')
		files_dict[(i, j)] = fout1
	for i in chr_list:
		fout = open(file_prefix + i + "." + i, 'w')
		files_dict[(i, i)] = fout
	return files_dict

def is_sorted_queryname(header):
	"""
	Check if bam fiel is sorted by read name.
	"""
	if("HD" in header):
		if("SO" in header["HD"]):
			if(header["HD"]["SO"] == "queryname"):
				return True
	return False

def extract_chr_list(header):
	chrs = []
	for chr_line in header["SQ"]:
		if "SN" in chr_line:
			chrs.append(chr_line["SN"])
	return(chrs)

main()
