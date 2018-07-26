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
import pybedtools
import time
import os
import subprocess
import glob
import logger

def split_main(input_bam, outdir, prefix, cutoff, per_chr, generate_hic):
	#sys.stdout = logger.Logger(outdir + "/" + prefix + ".feather.log")
	print(time.ctime() + " starting the splitting operation")
	samfile = pysam.AlignmentFile(input_bam, "rb")
	chr_list = extract_chr_list(samfile.header)
	autosomal_chrs = [chr_name for chr_name in chr_list if (chr_name.find('_') == -1 and chr_name[3].isdigit())]
	print(autosomal_chrs)
	if (len(chr_list) > 100):
		print("NOTE: Too many chromosomes. Generating bed files only for ones without an underscore in their name")
		chr_list = [chr_name for chr_name in chr_list if (chr_name.find('_') == -1)]
	if not is_sorted_queryname(samfile.header):
		sys.exit("Error: bam needs to be sorted by read name")
	filename_prefix = outdir + "/" + prefix
	temp_filename_prefix = outdir + "/tempfiles/" + prefix
	fname_output_shrt = temp_filename_prefix + ".shrt.bam"
	if (per_chr):
		short_bed_files, long_intra_bam_files, long_intra_bam_filenames, long_intra_bedpe_files = generate_per_chrom_files(chr_list, 
									filename_prefix, temp_filename_prefix, samfile.header)
	fname_output_long_intra = temp_filename_prefix + ".long.intra.bam"
	fout_long_intra = pysam.AlignmentFile(fname_output_long_intra, "wb", header=samfile.header)
	fname_long_intra_bedpe = filename_prefix + ".long.intra.bedpe"
	fname_output_long_inter = temp_filename_prefix + ".long.inter.bam"
	fname_output_shrt_inter = temp_filename_prefix + ".shrt.inter.bam"
	fname_shrt_bed = filename_prefix + ".shrt.vip.bed"
	fout_long_inter = pysam.AlignmentFile(fname_output_long_inter, "wb", header=samfile.header)
	fout_shrt_inter = pysam.AlignmentFile(fname_output_shrt_inter, "wb", header=samfile.header)
	fout_shrt = pysam.AlignmentFile(fname_output_shrt, "wb", header=samfile.header)
	fout_shrt_bed = open(fname_shrt_bed, "w")
	hic_chr_list = chr_list
	if (len(hic_chr_list) > 30):
		print("NOTE: Too many chromosomes. Generating hic only for ones without an underscore in their name")
		hic_chr_list = [chr_name for chr_name in hic_chr_list if (chr_name.find('_') == -1)]
	if (generate_hic):
		hic_tsv_files = open_hic_files(hic_chr_list, outdir, prefix)
	read_num = 0
	prev = pysam.AlignedSegment()
	prev.query_name = ""
	shrt_count = 0
	short_inter_count = 0
	shrt_vip_count = 0
	long_intra_count = 0
	long_inter_count = 0
	long_autosomal_intra_count = 0
	long_filtered_count = 0
	short_autosomal_count = 0
	short_auto_vip_count = 0
	short_autosomal_inter_count = 0
	for read in samfile:
		if read_num % 2 == 1:
			if read.query_name != prev.query_name:
				exit("unmatched query" + read.query_name)
		pos1 = prev.reference_end if prev.is_reverse else prev.reference_start
		pos2 = read.reference_end if read.is_reverse else read.reference_start

		if ((read.next_reference_name == read.reference_name) and (abs(read.template_length) <= cutoff) and (abs(read.template_length) >= 0)):
			fout_shrt.write(read)
			shrt_count += 1
			if (read.reference_name in autosomal_chrs):
				short_autosomal_count += 1
			if read_num % 2 == 1:
				if((prev.is_reverse == False and read.is_reverse == True) or (prev.is_reverse == True and read.is_reverse == False)):
					chrom = prev.reference_name
					if(pos1 != pos2):
						fout_shrt_bed.write("\t".join([chrom] + list(map(str, sorted([pos1, pos2]))) + [prev.query_name, "\n"]))
						shrt_vip_count += 1
						if chrom in autosomal_chrs:
							short_auto_vip_count += 1
						if (per_chr and (chrom in chr_list)):
							short_bed_files[chrom].write("\t".join([chrom] + list(map(str, sorted([pos1, pos2]))) + [prev.query_name, "\n"]))
		elif (read.next_reference_name == read.reference_name):
			if (per_chr and (read.reference_name in chr_list)):
				long_intra_bam_files[read.reference_name].write(read)
			if (read.reference_name in autosomal_chrs):
				long_autosomal_intra_count += 1
				if (abs(read.template_length) <= 1000000):
					long_filtered_count += 1
			fout_long_intra.write(read)
			long_intra_count += 1
		elif (abs(read.template_length) <= cutoff) and (abs(read.template_length) >= 0):
			if (read.reference_name in autosomal_chrs and prev.reference_name in autosomal_chrs):
				short_autosomal_inter_count += 1
			fout_shrt_inter.write(read)
			short_inter_count += 1
		else:
			long_inter_count += 1
			fout_long_inter.write(read)
		if read_num % 2 == 1:
			if (generate_hic):
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
				if (chr1 in hic_chr_list) and (chr2 in hic_chr_list):
					hic_tsv_files[tuple(sorted([chr1, chr2]))].write("\t".join([r1.query_name, str(strand1), chr1, str(r1.reference_start), "0", str(strand2), chr2, str(r2.reference_start), "1", "60", "60\n"]))
		read_num += 1
		prev = read
	fout_shrt_bed.close()
	samfile.close()
	fout_long_inter.close()
	fout_shrt_inter.close()
	fout_shrt.close()
	if (per_chr):
		close_files(short_bed_files)
		close_files(long_intra_bam_files)
	if (generate_hic):
		close_files(hic_tsv_files)
		hic_output =  outdir + "/" + prefix + ".hic.input"
		temp_files = glob.glob( outdir + "/tempfiles/hic_tempfiles/" + prefix + ".hic.*")
		command = ["cat"]
		command.extend(temp_files)
		#print(command)
		with open(hic_output, 'w') as outfile:
			
			proc = subprocess.Popen(command, stdout = outfile, shell= False)
			proc.communicate()
#		for chrom in chr_list:
#			short_bed_files[chrom].close()
#			long_intra_bam_files[chrom].close()
	fout_long_intra.close()
	print(time.ctime() + " generating long, intra-chromosomal bedpe file(s)")
	if (per_chr):
		for chrom in chr_list:
			sam = pysam.AlignmentFile(long_intra_bam_filenames[chrom])
			if (any(True for _ in sam)):
				long_bam = pybedtools.BedTool(long_intra_bam_filenames[chrom])
				long_bedpe = long_bam.bam_to_bed(bedpe = True, stream = True)
				long_bedpe.saveas(long_intra_bedpe_files[chrom])
			else:
				print("\tno bedpe generated for" + chrom +". Empty bam file.")
	print(time.ctime() + " writing to the combined bedpe file")
	long_bam = pybedtools.BedTool(fname_output_long_intra)
	long_bedpe = long_bam.bam_to_bed(bedpe = True, stream = True)
	long_bedpe.saveas(fname_long_intra_bedpe)
	qc_filename = outdir + "/" + prefix + ".feather.qc"
	read_count = read_num / 2
	long_intra_count /= 2
	shrt_count /= 2
	short_autosomal_count /= 2
	long_autosomal_intra_count /= 2
	short_autosomal_inter_count /= 2
	long_filtered_count /= 2
	short_inter_count /= 2
	long_inter_count /= 2
	inter_all_count = short_inter_count + long_inter_count
	with open(qc_filename, 'a') as outfile:
		outfile.write("{0:70} {1} ".format("number of intrachromosomal pairs", str(int(shrt_count + long_intra_count))))
		outfile.write("\t({0:.2f}%)\n".format(100 * int(float(shrt_count + long_intra_count)) / int(float(read_count))))
		outfile.write("{0:70} {1} ".format("number of short-range intrachromosomal pairs", str(int(shrt_count))))
		outfile.write("\t({0:.2f}%)\n".format(100 * int(float(shrt_count)) / int(float(read_count))))
		outfile.write("{0:70} {1} ".format("number of short-range intrachromosomal autosomal pairs", str(int(short_autosomal_count))))
		outfile.write("\t({0:.2f}%)\n".format(100 * int(float(short_autosomal_count)) / int(float(read_count))))
		outfile.write("{0:70} {1} ".format("number of short-range autosomal vip pairs", str(int(short_auto_vip_count))))
		outfile.write("\t({0:.2f}%)\n".format(100 * int(float(short_auto_vip_count)) / int(float(read_count))))
		outfile.write("{0:70} {1} ".format("number of long-range intra-chromosomal pairs:", str(int(long_intra_count))))
		outfile.write("\t({0:.2f}%)\n".format(100 * int(float(long_intra_count)) / int(float(read_count))))
		outfile.write("{0:70} {1} ".format("number of inter-chromosomal pairs:", str(int(inter_all_count))))
		outfile.write("\t({0:.2f}%)\n".format(100 * int(float(inter_all_count)) / int(float(read_count))))
		outfile.write("{0:70} {1} ".format("number of long-range autosomal intra-chromosomal pairs", str(int(long_autosomal_intra_count))))
		outfile.write("\t({0:.2f}%)\n".format(100 * int(float(long_autosomal_intra_count)) / int(float(read_count))))
		outfile.write("{0:70} {1} ".format("number of long-range autosomal intra 1k-1M pairs", str(int(long_filtered_count))))
		outfile.write("\t({0:.2f}%)\n".format(100 * int(float(long_filtered_count)) / int(float(read_count))))
	print(time.ctime() + " splitting completed")

def close_files(files_dict):
	for filename in files_dict.values():
		filename.close()

def open_hic_files(chr_list, outdir, prefix):
	files_dir =  outdir + "/tempfiles/hic_tempfiles"
	if not os.path.exists(files_dir):
		os.mkdir(files_dir, 0o777)
	file_prefix = files_dir + "/" + prefix + ".hic."
	files_dict = {}
	for i, j in combinations(chr_list, 2):
		i, j = sorted([i, j])
		fout1 = open(file_prefix + i + "." + j, 'w')
		files_dict[(i, j)] = fout1
	for i in chr_list:
		fout = open(file_prefix + i + "." + i, 'w')
		files_dict[(i, i)] = fout
	return files_dict

def generate_per_chrom_files(chr_list, filename_prefix, temp_filename_prefix, header):
	bam_files = {}
	bam_filenames = {}
	bedpe_files = {}
	short_bed_files = {}
	for chrom in chr_list:
		fname_output_long_intra = temp_filename_prefix + "." + chrom + ".long.intra.bam"
		fname_output_short_bed = filename_prefix + "." + chrom + ".shrt.vip.bed"
		fout_short_bed = open(fname_output_short_bed, "w")
		fout_bam = pysam.AlignmentFile(fname_output_long_intra, "wb", header = header)
		fname_bedpe = filename_prefix + "." + chrom + ".long.intra.bedpe"
		short_bed_files[chrom] = fout_short_bed
		bam_filenames[chrom] = fname_output_long_intra
		bam_files[chrom] = fout_bam
		bedpe_files[chrom] = fname_bedpe
	return (short_bed_files, bam_files, bam_filenames, bedpe_files)

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

