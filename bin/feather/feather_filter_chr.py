import subprocess
import argparse
import os.path
import pandas
import pysam
from io import StringIO
import itertools
import numpy as np
import random
import re
#import pysam
from tempfile import TemporaryFile
import time
import copy
import logger
import sys
from shutil import copyfile

sam_columns = {"query_name" : 0, "flag" : 1, "chr_name" : 2, "pos" : 3, "mapq" : 4,
		"cigar" : 5, "mate_name" : 6, "mate_pos" : 7, "tlen" : 8, "seq" : 9, "quality" : 10}

flag_table_proper={(0, 16): (99, 147), (2048, 16): (99, 147), (0, 2064):  (99, 147), (2048, 2064): (99, 147),     # + -
			(16, 0): (83, 163), (16, 2048): (83, 163), (2064, 0):  (83, 163), (2064, 2048): (83, 163),     # - +
			(0,  0): (65, 129), (0,  2048): (65, 129), (2048, 0):  (65, 129), (2048, 2048): (65, 129),     # + +
			(16, 16):(113, 177), (16, 2064):(113, 177), (2064, 16):(113, 177), (2064, 2064):(113, 177),    # - -
		}
#qseq="K00168:100:HF52YBBXX:4:1101:27580:1103"
qseq = "K00168:100:HF52YBBXX:4:1101:27336:1103"

def filter_main(fastq1, fastq2, bwa_index, mapq, outdir, prefix, threads, to_file = False):
	sys.stdout = logger.Logger(outdir + "/" + prefix + ".feather.log")
	print(time.ctime() + " starting mapping and filtering operation")
	check_arguments(fastq1, fastq2, bwa_index, mapq, threads)
	paired_filename, bwa1_filename, bwa2_filename, bwa1_sorted_filename, bwa2_sorted_filename, combined_bwa_filename, qc_filename = set_filenames(fastq1, fastq2, outdir, prefix)
	#running bwa mem
	for fastq, bwa_filename in [(fastq1, bwa1_filename), (fastq2, bwa2_filename)]:
		if fastq.endswith(".fastq") or fastq.endswith("fastq.gz") or fastq.endswith("fq") or fastq.endswith("fq.gz"):
			bwa_mem(fastq, bwa_index, threads, bwa_filename)
		elif not (fastq.endswith(".sam") or fastq.endswith(".bam")):
			exit("Error: Input file for filtering should be of type fastq, fastq.gz, sam, or bam. Exiting!")
	if bwa1_filename.endswith(".bam"):
		proc = subprocess.Popen("samtools view " + bwa1_filename + " | awk ' $1 !~ /@/ {print $1}' " + "| uniq -c|wc -l", stdout = subprocess.PIPE, shell = True)
		read_count = proc.stdout.read().decode("utf-8")
	else:
		proc = subprocess.Popen("awk ' $1 !~ /@/ {print $1}' " + bwa1_filename  + "| uniq -c|wc -l", stdout = subprocess.PIPE, shell = True)
		read_count = proc.stdout.read().decode("utf-8")

	#pairing and filtering alignments for chimeric reads
	for bwa_filename, bwa_sorted_filename in ([bwa1_filename, bwa1_sorted_filename], [bwa2_filename, bwa2_sorted_filename]):
		bwa = pysam.AlignmentFile(bwa_filename)
		if not is_sorted_queryname(bwa.header):
			print(time.ctime() + " calling samtools sort for " + bwa_filename + " storing in " + bwa_sorted_filename)
			pysam.sort("-o", bwa_sorted_filename , "-n", "-@", str(threads), bwa_filename)
		else:
			copyfile(bwa_filename, bwa_sorted_filename)
	print(time.ctime() + " merging " + bwa1_sorted_filename + " and " + bwa2_sorted_filename)
	pysam.merge("-n", "-f",  combined_bwa_filename, bwa1_sorted_filename, bwa2_sorted_filename)
	print(time.ctime() + " filtering and pairing reads")
	filter_pair_reads(combined_bwa_filename, mapq, paired_filename, qc_filename)
	print(time.ctime() + " paired bam file generated. Sorting by coordinates.")
	pysam.sort("-o", paired_filename + ".srt.bam", "-@", str(threads), paired_filename + ".bam")
	print(time.ctime() + " calling samtools rmdup")
	pysam.rmdup(paired_filename + ".srt.bam", paired_filename + ".rmdup.bam")
	#proc = subprocess.Popen(["samtools", "rmdup", paired_filename + ".srt.bam", paired_filename + ".rmdup.bam"])
	#proc.communicate()
	print(time.ctime() + " calling samtools flagstat on mapped file")
	proc = subprocess.Popen("samtools flagstat " + paired_filename + ".srt.bam > " + paired_filename + ".srt.bam.flagstat", 
				shell = True)
	proc.communicate()
	with open(paired_filename + ".srt.bam.flagstat") as flag_file:
		lines = flag_file.readlines()
		uniquely_mapped_count = lines[7].split()[0]
	print(time.ctime() + " calling samtools flagstat on mapped and duplicate-removed file")
	proc = subprocess.Popen("samtools flagstat " + paired_filename + ".rmdup.bam > " + paired_filename + ".rmdup.flagstat", 
				shell = True)
	proc.communicate()
	with open(paired_filename + ".rmdup.flagstat") as flag_file:
		lines = flag_file.readlines()
		duprmd_count = lines[7].split()[0]
		intra_count = lines[11].split()[0]
		intra_count = str(int(float(intra_count)) / 2)
	print(time.ctime() + " calling samtools sort for sorting by query names")
	#pysam.sort("-n", "-o", bwa_filename + ".srtn.rmdup.bam", paired_filename + ".rmdup.bam")
	pysam.sort("-o", paired_filename + ".srtn.rmdup.bam", 
				"-@", str(threads), "-n", paired_filename + ".rmdup.bam")
	#proc.communicate()
	#proc.wait()
	print(time.ctime() + " finishing filtering")
	qc_filename = outdir + "/" + prefix + ".feather.qc"
	with open(qc_filename, 'w') as outfile:
		outfile.write("{0:70} {1}".format("number of sequencing pairs", str(read_count)))
		outfile.write("{0:70} {1} ".format("number of uniquely mapped pairs (MAPQ >= " + str(mapq) + ")", str(uniquely_mapped_count)))
		outfile.write("\t({0:.2f}%)\n".format(100 * (int(float(uniquely_mapped_count)) / int(float(read_count)))))
		outfile.write("{0:70} {1} ".format("number of pairs after duplicate removal", str(duprmd_count)))
		outfile.write("\t({0:.2f}%)\n".format(100 * (int(float(duprmd_count)) / int(float(read_count)))))
		#outfile.write("{0:70} {1} ".format("number of interchromosomal pairs", str(intra_count)))
		#outfile.write("\t({0:.2f}%)\n".format(100 * int(float(intra_count)) / int(float(read_count))))
	return (paired_filename + ".srtn.rmdup.bam")

def is_sorted_queryname(header):
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

def filter_pair_reads(bwa_filename, mapq, paired_filename, qc_filename):
	merged_sam = pysam.AlignmentFile(bwa_filename, "rb")
	with pysam.AlignmentFile(paired_filename + ".bam", "wb", header = merged_sam.header) as outfile:
		current = ""
		counter = 0
		aligned_reads = []
		paired_reads = []
		unmapped_read_count = 0
		total_read_count = 0
		chimeric_read_count = 0
		chrs = set()
		for read in merged_sam.fetch(until_eof = True):
			total_read_count += 1
			if read.flag == 4 or read.mapq < mapq:
				unmapped_read_count += 1
				continue
			else:
				if read.query_name != current:
					if current != "":
						if counter == 2:
							paired_reads.append(pair_up(aligned_reads))
						elif counter == 3:
							chimeric_read_count += 1
							chr_count = len(chrs)
							if chr_count == 1 or chr_count == 2:
								selected_pair = select_valid_pair(aligned_reads, counter, chr_count)
								paired_reads.append(pair_up([aligned_reads[i] for i in selected_pair]))
						for i, j in paired_reads:
							outfile.write(i)
							outfile.write(j)
						paired_reads = []
					aligned_reads = [read]
					counter = 1
					chrs = set()
					chrs.add(read.reference_name)
					current = read.query_name
				else:
					counter += 1
					chrs.add(read.reference_name)
					aligned_reads.append(read)
		if counter == 2:
			paired_reads.append(pair_up(aligned_reads))
		elif counter == 3:
			chimeric_read_count += 1
			selected_pair = select_valid_pair(aligned_reads, counter, len(chrs))
			paired_reads.append(pair_up([aligned_reads[i] for i in selected_pair]))
		for i, j in paired_reads:
			outfile.write(i)
			outfile.write(j)
	with open(qc_filename, "a") as qc_fout:
		output = " ".join(["Total reads:", str(total_read_count)])
		output += "\n" + " ".join(["Unmapped/low quality reads:", str(unmapped_read_count)])
		output += "\n" + " ".join(["Chimeric reads (after low quality read removal):", str(chimeric_read_count)])
		qc_fout.write(output)

def set_filenames(fastq1, fastq2, outdir, prefix):
	tempdir = outdir + "/tempfiles"
	if not os.path.exists(tempdir):
		os.makedirs(tempdir)
	paired_filename = outdir + "/" + prefix + ".paired"
	if fastq1.endswith(".fastq") or fastq1.endswith(".fastq.gz") or fastq1.endswith(".fq") or fastq1.endswith(".fq.gz"):
		fastq1_prefix = fastq1[ (fastq1.rfind("/") + 1) : ]
		bwa1_filename = tempdir + "/" + fastq1_prefix + ".bwa.sam"
		bwa1_sorted_filename = bwa1_filename + ".srtn"
	elif fastq1.endswith(".sam") or fastq1.endswith(".bam"):
		setname = fastq1[ (fastq1.rfind("/") + 1) : ]
		bwa1_filename = fastq1
		bwa1_sorted_filename = tempdir + "/" + setname + ".srtn"
	if fastq2.endswith(".fastq") or fastq2.endswith(".fastq.gz") or fastq1.endswith(".fq") or fastq1.endswith(".fq.gz"):
		fastq2_prefix = fastq2[ (fastq2.rfind("/") + 1) : ]
		bwa2_filename = tempdir + "/" + fastq2_prefix + ".bwa.sam"
		bwa2_sorted_filename = bwa2_filename + ".srtn"
	elif fastq2.endswith(".sam") or fastq2.endswith(".bam"):
		setname = fastq2[ (fastq2.rfind("/") + 1) : ]
		bwa2_filename = fastq2
		bwa2_sorted_filename = tempdir + "/" + setname + ".srtn"
	combined_bwa_filename = tempdir + "/" + prefix + ".merged.srtn.bam"
	qc_filename = outdir + "/" + prefix + ".qc"
	return(paired_filename, bwa1_filename, bwa2_filename, bwa1_sorted_filename, bwa2_sorted_filename, combined_bwa_filename, qc_filename)

def set_tempfile(input_content = None, output_content = None, binary = True):
	tfile = TemporaryFile("wb") if (binary) else TemporaryFile("w")
	if (input_content):
		tfile.write(input_content)
		tfile.seek(0)
	return(tfile)
		

def bwa_mem(fastq, bwa_index, threads, output_filename):
	print(time.ctime() + " calling bwa for " + fastq)
	output_file = open(output_filename, "w")
	proc = subprocess.Popen(["bwa", "mem", "-t", str(threads), bwa_index, fastq], stdout = output_file, stderr = open(output_filename + ".log", 'w'))
	proc.wait()
	output_file.close()

def select_valid_pair(aligned_reads, count, chr_count):
	if (chr_count == 1):
		combinations = [(0, 1), (0, 2), (1, 2)]
		dists = []
		for i, j in combinations:
			dists.append(abs(aligned_reads[i].reference_start - aligned_reads[j].reference_start))
		arg_order = np.argsort(dists)
		second_largest_pair = combinations[arg_order[1]]
		pair_indices = second_largest_pair
	elif (chr_count == 2):
		if aligned_reads[0].reference_name == aligned_reads[1].reference_name:
			base_index = 2
			matched_indices = [0, 1]
		elif aligned_reads[0].reference_name == aligned_reads[2].reference_name:
			base_index = 1
			matched_indices = [0, 2]
		else:
			base_index = 0
			matched_indices = [1, 2]
		pair_indices = sorted([base_index, min(matched_indices)])
	return (pair_indices)
		
def pair_up(read_pair):
	r1_cp = pysam.AlignedSegment()
	r2_cp = pysam.AlignedSegment()
	r1_cp = copy.deepcopy(read_pair[0])
	r2_cp = copy.deepcopy(read_pair[1])
	if r1_cp.query_name != r2_cp.query_name:
		print("Error: read name unmathced.\n")
		sys.exit(1);
	# change flag
	flag_swag = flag_table_proper[(r1_cp.flag, r2_cp.flag)]
	r1_cp.flag = flag_swag[0]
	r2_cp.flag = flag_swag[1]
	# now change RNEXT and PNEXT
	if(r1_cp.reference_name  == r2_cp.reference_name):
		#r1_cp.next_reference_name = r1_cp.reference_name
		#r2_cp.next_reference_name = r1_cp.reference_name
		r1_cp.next_reference_start = r2_cp.reference_start
		r2_cp.next_reference_start = r1_cp.reference_start
		r1_cp.template_length  = r1_cp.next_reference_start -  r1_cp.reference_start
		r2_cp.template_length  = -r1_cp.template_length
		r1_cp.next_reference_name = r2_cp.next_reference_name = "="
	else:
		r1_cp.next_reference_name, r2_cp.next_reference_name = r2_cp.reference_name, r1_cp.reference_name
		r1_cp.next_reference_start = r2_cp.reference_start
		r2_cp.next_reference_start = r1_cp.reference_start
		r1_cp.template_length  = r1_cp.next_reference_start -  r1_cp.reference_start
		r2_cp.template_length  = -r1_cp.template_length
	return (r1_cp, r2_cp)

def check_requirements():
	for requirement in requires:
		status_which, result = subprocess.getstatusoutput("which " + requirement)
		status_ls = os.path.exists("utils/" + requirement)
		if (status_which % 256 != 0 and not status_ls):
			exit("ERROR: " + requirement + " cannot be found. Please make sure it is " +
			"installed and is in the system path. Exiting!")

def check_arguments(fastq1, fastq2, bwa_index, mapq, threads):
	check_file_existance((fastq1, fastq2))
	if fastq1 == fastq2: exit("ERROR: Input fastq files should be different. Exiting!")
	if (int(threads) < 1):
		print(threads)
		exit("ERROR: Number of threads (-t) should be a positive integer. Exiting!")

def check_file_existance(files):
	for file_name in files:
		if not os.path.exists(file_name):
			exit("ERROR: Input file " + file_name + " does not exist. Please make sure " +
			"the correct path is provided. Exiting!")

