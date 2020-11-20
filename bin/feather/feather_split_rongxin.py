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

def split_main(input_bam, outdir, prefix, cutoff, per_chr, generate_hic, chip_peaks):
	#sys.stdout = logger.Logger(outdir + "/" + prefix + ".feather.log")
	print(time.ctime() + " starting the splitting operation")
	flagstat_filename = glob.glob(outdir + "/*.rmdup.flagstat")
	if chip_peaks:
		chip_peaks = glob.glob(chip_peaks)
	else:
		chip_peaks = ""
	if len(flagstat_filename) > 0 and len(chip_peaks) > 0:
		flagstat_filename = flagstat_filename[0]
		chip_peaks = chip_peaks[0]
		print(time.ctime() + " QC metrics will be computed from: " + chip_peaks + " and " + flagstat_filename)
	elif len(flagstat_filename) > 0:
		chip_peaks = False
		flagstat_filename = flagstat_filename[0]
		print(time.ctime() + " WARNING: no ChIP peaks file found to compute the on-chip ratio")
		print(time.ctime() + " QC metrics will be computed from: " + flagstat_filename)
	elif len(chip_peaks) > 0:
		if len(flagstat_filename) == 0:
			flagstat_filename = input_bam + ".rmdup.flagstat"
			print(time.ctime() + " calling samtools flagstat on provided bam file")
			proc = subprocess.Popen("samtools flagstat " + input_bam + " > " + flagstat_filename, shell = True)
			proc.communicate()
		chip_peaks = chip_peaks[0]
		print(time.ctime() + " QC metric (on-chip ratio) will be computed from: " + chip_peaks + " and " + flagstat_filename)
	samfile = pysam.AlignmentFile(input_bam, "rb")
	filename_prefix = outdir + "/" + prefix
	temp_filename_prefix = outdir + "/tempfiles/" + prefix
	chr_list = extract_chr_list(samfile.header)
	chrname_lengths = list(map(len, chr_list))
	rename_chrs = True if (min(chrname_lengths) == 1) else False
	if rename_chrs:
		chr_list =["chr" + chrom for chrom in chr_list]
		s2 = pysam.AlignmentFile(temp_filename_prefix + "_sam_header", "wb", reference_names = chr_list, reference_lengths = samfile.header.lengths)
		new_header = s2.header
	else:
		new_header = samfile.header
	autosomal_chrs = [chr_name for chr_name in chr_list if ((chr_name.find('_') == -1 and chr_name.find('.') == -1 and chr_name[3].isdigit()) or chr_name == "chrX" or chr_name == "chrY")]
	potential_chrs = ['chr' + str(x) for x in list(range(1,23)) + ['X', 'Y']]
	autosomal_chrs = [x for x in autosomal_chrs if x in potential_chrs]
	chr_list = autosomal_chrs
	if not is_sorted_queryname(samfile.header):
		sys.exit("Error: bam needs to be sorted by read name")
	filename_prefix = outdir + "/" + prefix
	temp_filename_prefix = outdir + "/tempfiles/" + prefix
	fname_output_shrt = temp_filename_prefix + ".shrt.bam"
	#fname_output_shrt_genomewide =  temp_filename_prefix + ".shrt.genomewide.bam"
	if (per_chr):
		short_bed_files, long_intra_bam_files, long_intra_bam_filenames, long_intra_bedpe_files = generate_per_chrom_files(chr_list, 
									filename_prefix, temp_filename_prefix, new_header)
	fname_output_long_intra = temp_filename_prefix + ".long.intra.bam"
	fout_long_intra = pysam.AlignmentFile(fname_output_long_intra, "wb", header = new_header)
	fname_long_intra_bedpe = filename_prefix + ".long.intra.bedpe"
	fname_output_long_inter = temp_filename_prefix + ".long.inter.bam"
	fname_output_shrt_inter = temp_filename_prefix + ".shrt.inter.bam"
	fname_shrt_bed = filename_prefix + ".shrt.vip.bed"
	fout_long_inter = pysam.AlignmentFile(fname_output_long_inter, "wb", header = new_header)
	fout_shrt_inter = pysam.AlignmentFile(fname_output_shrt_inter, "wb", header = new_header)
	fout_shrt = pysam.AlignmentFile(fname_output_shrt, "wb", header = new_header)
	fout_shrt_bed = open(fname_shrt_bed, "w")
	hic_chr_list = autosomal_chrs
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
		if rename_chrs:
			read = rename_read_chromosomes(read, new_header)
			#autosomal_chrs = ["chr" + chrom for chrom in autosomal_chrs]
			#print(autosomal_chrs)
		if read_num % 2 == 1:
			if read.query_name != prev.query_name:
				exit("unmatched query" + read.query_name)
		pos1 = prev.reference_end if prev.is_reverse else prev.reference_start
		pos2 = read.reference_end if read.is_reverse else read.reference_start

		is_intra, is_short, is_vip, is_auto = classify_reads(read, prev, autosomal_chrs, cutoff)
		if is_intra and is_auto and (not is_short): #intra-chromosomal long autosomal
			long_autosomal_intra_count += 1
			if (per_chr):
				long_intra_bam_files[read.reference_name].write(read)
			if (abs(read.template_length) <= 1000000):
				long_filtered_count += 1
			fout_long_intra.write(read)
			long_intra_count += 1
		if (not is_intra): #inter-chromosomal autosomal reads
			long_inter_count += 1
			if is_auto:
				fout_long_inter.write(read)
		if (is_short and is_intra): #short reads
			if is_auto: # write to bam file
				fout_shrt.write(read)
				shrt_count += 1
			if is_vip and is_auto and pos1 != pos2: #short vip reads
				if read_num % 2 == 1:
					chrom = prev.reference_name
					fout_shrt_bed.write("\t".join([chrom] + list(map(str, sorted([pos1, pos2]))) + [prev.query_name, "\n"]))
					shrt_vip_count += 1
					if (per_chr):
						short_bed_files[chrom].write("\t".join([chrom] + list(map(str, sorted([pos1, pos2]))) + [prev.query_name, "\n"]))
		if read_num % 2 == 1:
			if (generate_hic):
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
	#inter_all_count = short_inter_count + long_inter_count
	intra_all_count = shrt_count + long_intra_count
	if (flagstat_filename):
		with open(flagstat_filename) as flag_file:
			lines = flag_file.readlines()
			duprmd_count = int(lines[7].split()[0])
			trans_ratio = long_inter_count / duprmd_count
			long_cis_ratio = long_intra_count / intra_all_count
	else:
		trans_ratio, long_cis_ratio = -1, -1
	if (chip_peaks):
		command = "cut -f1-3 " + fname_shrt_bed + " > " + outdir + "/tempfiles/shrt_vip_bed.cut"
		proc = subprocess.Popen(command, stdout = subprocess.PIPE, shell= True)
		proc.communicate()
		command = "cut -f1-3 " + chip_peaks + " > " + outdir + "/tempfiles/chip_peaks.cut"
		proc = subprocess.Popen(command, stdout = subprocess.PIPE, shell= True)
		proc.communicate()
		command = "bedtools intersect -a " + outdir + "/tempfiles/shrt_vip_bed.cut" + " -b " + outdir + "/tempfiles/chip_peaks.cut" + " -u | wc -l"
		proc = subprocess.Popen(command, stdout = subprocess.PIPE, shell= True)
		on_chip = int(proc.stdout.read().decode("utf-8"))
		command = "wc -l " + outdir + "/tempfiles/shrt_vip_bed.cut"
		proc = subprocess.Popen(command, shell = True, stdout = subprocess.PIPE)
		chip_count = int(proc.stdout.read().decode("utf-8").split()[0])
		on_chip_ratio = float(on_chip) / chip_count
	else:
		on_chip_ratio = -1
	with open(qc_filename, 'a') as outfile:
		outfile.write("{0:70} {1} ".format("number of intrachromosomal pairs", str(int(intra_all_count))))
		outfile.write("\t({0:.2f}%)\n".format(100 * int(float(intra_all_count)) / int(float(read_count))))
		outfile.write("{0:70} {1} ".format("number of short-range intrachromosomal pairs", str(int(shrt_count))))
		outfile.write("\t({0:.2f}%)\n".format(100 * int(float(shrt_count)) / int(float(read_count))))
		#outfile.write("{0:70} {1} ".format("number of short-range intrachromosomal autosomal pairs", str(int(short_autosomal_count))))
		#outfile.write("\t({0:.2f}%)\n".format(100 * int(float(short_autosomal_count)) / int(float(read_count))))
		outfile.write("{0:70} {1} ".format("number of short-range vip pairs", str(int(shrt_vip_count))))
		outfile.write("\t({0:.2f}%)\n".format(100 * int(float(shrt_vip_count)) / int(float(read_count))))
		outfile.write("{0:70} {1} ".format("number of long-range intrachromosomal pairs", str(int(long_intra_count))))
		outfile.write("\t({0:.2f}%)\n".format(100 * int(float(long_intra_count)) / int(float(read_count))))
		outfile.write("{0:70} {1} ".format("number of interchromosomal pairs", str(int(long_inter_count))))
		outfile.write("\t({0:.2f}%)\n".format(100 * int(float(long_inter_count)) / int(float(read_count))))
		#outfile.write("{0:70} {1} ".format("number of long-range autosomal intra-chromosomal pairs", str(int(long_autosomal_intra_count))))
		#outfile.write("\t({0:.2f}%)\n".format(100 * int(float(long_autosomal_intra_count)) / int(float(read_count))))
		#outfile.write("{0:70} {1} ".format("number of long-range autosomal intra 1k-1M pairs", str(int(long_filtered_count))))
		#outfile.write("\t({0:.2f}%)\n".format(100 * int(float(long_filtered_count)) / int(float(read_count))))
		outfile.write("{0:70}".format("trans read ratio") + "{0:.5f} ".format(trans_ratio))
		outfile.write("\t({0:.2f}%)\n".format(100 * trans_ratio))
		outfile.write("{0:70}".format("long cis read ratio") + "{0:.5f} ".format(long_cis_ratio))
		outfile.write("\t({0:.2f}%)\n".format(100 * long_cis_ratio))
		outfile.write("{0:70}".format("FRiP") + "{0:.5f} ".format(on_chip_ratio))
		outfile.write("\t({0:.2f}%)\n".format(100 * on_chip_ratio))
	print(time.ctime() + " splitting completed")


def classify_reads(read1, read2, autosomal_chrs, cutoff):
	is_intra, is_short, is_vip, is_auto = [False] * 4
	if (read1.reference_name in autosomal_chrs) and (read1.next_reference_name in autosomal_chrs):
		is_auto = True
	if ((abs(read1.template_length) <= cutoff) and (abs(read1.template_length) >= 0)):
		is_short = True
	if (read1.next_reference_name == read1.reference_name):
		is_intra = True
	if ((read2.is_reverse == False and read1.is_reverse == True) or (read2.is_reverse == True and read1.is_reverse == False)):
		is_vip = True
	return is_intra, is_short, is_vip, is_auto

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

def rename_read_chromosomes(read, header):
	r = pysam.AlignedSegment(header)
	r.reference_name = "chr" + read.reference_name
	r.next_reference_name = "chr" + read.next_reference_name
	r.tags = read.tags
	r.bin, r.cigar, r.cigarstring, r.cigartuples, r.flag, r.is_duplicate, r.is_paired, r.is_proper_pair, r.is_qcfail, r.is_read1, r.is_read2, r.is_reverse, r.is_secondary, r.is_supplementary, r.is_unmapped, r.isize, r.mapping_quality, r.mapq, r.mate_is_reverse, r.mate_is_unmapped, r.mpos, r.mrnm, r.next_reference_id, r.next_reference_start, r.pnext, r.pos = read.bin, read.cigar, read.cigarstring, read.cigartuples, read.flag, read.is_duplicate, read.is_paired, read.is_proper_pair, read.is_qcfail, read.is_read1, read.is_read2, read.is_reverse, read.is_secondary, read.is_supplementary, read.is_unmapped, read.isize, read.mapping_quality, read.mapq, read.mate_is_reverse, read.mate_is_unmapped, read.mpos, read.mrnm, read.next_reference_id, read.next_reference_start, read.pnext, read.pos
	r.qname, r.query_name, r.query_sequence, r.reference_id, r.reference_start = read.qname, read.query_name, read.query_sequence, read.reference_id, read.reference_start
	r.rname, r.rnext, r.seq, r.tags, r.template_length, r.tid, r.tlen =  read.rname, read.rnext, read.seq, read.tags, read.template_length, read.tid, read.tlen
	r.query_qualities = read.query_qualities
	return(r)

if __name__ == "__main__":
	input_bam, outdir, prefix, cutoff, per_chr, generate_hic, chip_peaks = sys.argv[1:len(sys.argv)]
	split_main(input_bam, outdir, prefix, int(float(cutoff)), per_chr, generate_hic, chip_peaks)
