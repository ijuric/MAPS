#!/bin/bash
############################################################################
###                   Arima MAPS v1.8 QC Metrics Only                    ###
############################################################################

# Use this config file to run Arima HiChIP data, created using the Arima Hi-C+ kit, through the MAPS pipeline
# This pipeline was modified from MAPS: https://github.com/ijuric/MAPS
# All modifications to this code are in post processing steps to generate Arima specific quality control metrics and do not affect the MAPS algorithm itself.
# For technical assistance please contact Ming Hu at hum@ccf.org
# Please use the following citation for referencing MAPS: Juric I, Yu M, Abnousi A, Raviram R, Fang R, Zhao Y, et al. (2019) MAPS: Model-based analysis of long-range chromatin interactions from PLAC-seq and HiChIP experiments. PLoS Comput Biol 15(4): e1006982. https://doi.org/10.1371/ journal.pcbi.1006982

############################################################################
###                       Computing Environment                          ###
############################################################################

# Please install the following dependencies and add them tou your PATH variable to ensure MAPS will execute in your computing environment:

# Dependencies:

# Python 3.4 (or later): pandas (v0.20.3), numpy (v1.13.1), itertools (v3.2), pysam (v0.15.2), pybedtools (v0.8.0), deeptools
# R 3.4.3: argparse
# other: bedtools, samtools, and bcftools.

# see https://github.com/ijuric/MAPS/README.md for installation help

# please provide the path to the following tools:

python_path=[PYTHON_PATH]
Rscript_path=[Rscript_PATH]
java_path=[JAVA_PATH]
picard_path=[PICARD_PATH]
cwd=[MAPS_BIN_PATH]

############################################################################
###                          Customer Inputs                             ###
############################################################################

feather=1 # "1" to run feather, "0" to skip
maps=1 # "1" to run MAPS on data previously processed with feather, "0" to skip
# setting both feather and MAPS to zero will run the Arima QC metrics only on previous processed data
Arima_QC=1 # "1" to generate QC metrics, "0" to skip
number_of_datasets=1 # number of biological replicas. Can be 1,2, … N. If set to 1, script assumes that user is running MAPS for a particular biological replica. If number_of_datasets is larger than 1, script assumes that user is running MAPS on merged data, and it will look for .bedpe and .bed files of each replica defined in dataset1, dataset2,..., datasetN lines (See Additional Parameters section)
dataset_name=[DATASET_NAME] # sample name with no file extention, Ex: "sample1"
fastq_format=".fastq" # format that MAPS should expect the data in. Ex: "_001.fastq.gz"
pair_format="_" # fastq file names must be named [DATASET_NAME][pair_format][R1 or R2][fastq_format]. Ex: "_" or "."
fastq_dir=[FASTQ_DIR] # absolute path to the Directory the data is in. Ex: "/home/sample1/fastq/"
outdir=[OUTPUT_DIR] # directory which the feather and MAPS outputs will be placed
macs2_filepath=[CHIP_PEAK_FILE] # path to the ChIP peak file that will be used for loop calling with MAPS
organism=[REFERENCE_NAME] # organism of the genomic feature file to be used, options: "mm9", "mm10", "hg19", and "hg38"
bwa_index=[BWA_REFERENCE_SEQUENCE] # absolute path to the reference genome sequence (.fa) which will be used to derive the BWA indexes as well. Ex: "/home/reference_sequence/hg19.fa"
threads=[NUMBER_THREADS] # number of threads for MAPS to run on

# Specify the absolute path for genomics feature files below if not default
if [ $organism == "mm10" ]; then
    genomic_feat_filepath="$cwd/../Arima_Genomics/genomic_features/mm10_F_GC_M_arima_5Kb_el.txt"
    chr_count=19
elif [ $organism == "mm9" ]; then
    genomic_feat_filepath="$cwd/../Arima_Genomics/genomic_features/mm9_F_GC_M_arima_5Kb_el.txt"
    chr_count=19
elif [ $organism == "hg19" ]; then
    genomic_feat_filepath="$cwd/../Arima_Genomics/genomic_features/hg19_F_GC_M_arima_5Kb_el.txt"
    chr_count=22
elif [ $organism == "hg38" ]; then
    genomic_feat_filepath="$cwd/../Arima_Genomics/genomic_features/hg38_F_GC_M_arima_5Kb_el.txt"
    chr_count=22
fi

############################################################################
###                    Arima Recommended Parameters                      ###
############################################################################

# The parameter settings in this section are based on the default parameters for MAPS with some minor adjustments. These parameters have been optimized by internal benchmarking and have been found to optimize sensitivity while not inflating the false positive rate.

plot=1 # "1" to generate plots, "0" to skip
generate_hic=1 # "1" to generate a .hic file for visualization with Juicer, "0" to skip
bin_size=5000 # Resolution of the loops called.
binning_range=2000000 # Maximum distance for loop calling
fdr=2 # Do not change! False Discovery Rate threshold 1 = 0.1, 2 = 0.01, 3 = 0.001, ect...
filter_file="None"
mapq=30 # Phred scaled mapping quality threshold
length_cutoff=1000 # Minimum genomic distance
model="pospoisson" # Regression Model for Loop calling. Options are: "pospoisson" and "negbinom"
sex_chroms_to_process="NA" # Which Sex Chromosomes should be processed

############################################################################
###                       Additional Parameters                          ###
############################################################################

# SET THE VARIABLES AT THIS PORTION ONLY IF
# number_of_datasets > 1 (merging exisitng datasets)
# specify as many datasets as required

dataset1=""
dataset2=""
dataset3=""
dataset4=""
#...

###SET THESE VARIABLES ONLY IF FEATHER = 0 AND YOU WANT TO RUN
###USING A SPECIFIC FEATHER OUTPUT RATHER THAN $datasetname_Current

feather_output_symlink=""

############################################################################
###                           MAPS Pipeline                              ###
############################################################################

DATE=`date '+%Y%m%d_%H%M%S'`
fastq1=$fastq_dir/$dataset_name$pair_format"R1"$fastq_format
fastq2=$fastq_dir/$dataset_name$pair_format"R2"$fastq_format
feather_output=$outdir"/feather_output/"$dataset_name"_"$DATE"/"
if [ "$feather_output_symlink" == "" ]; then
	feather_output_symlink=$outdir"/feather_output/"$dataset_name"_current"
fi
per_chr='True' # set this to zero if you don't want per chromosome output bed and bedpe files
feather_logfile=$feather_output"/"$dataset_name".feather.log"
resolution=$(bc <<< "$bin_size/1000")
# cwd="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"
hic_dir="tempfiles/hic_tempfiles"

if [[ $sex_chroms_to_process != "X" && $sex_chroms_to_process != "Y" && $sex_chroms_to_process != "XY" ]]; then
	sex_chroms_to_processes="NA"
	sex_chroms=""
else
	sex_chroms=$sex_chroms_to_process
fi
long_bedpe_dir=$feather_output_symlink"/"
short_bed_dir=$feather_output_symlink"/"
maps_output=$outdir"/MAPS_output/"$dataset_name"_"$DATE"/"
maps_output_symlink=$outdir"/MAPS_output/"$dataset_name"_current"

if [ $feather -eq 1 ]; then
	mkdir -p $feather_output
	if [ $number_of_datasets -ge 2 ]; then
		ds_array=()
		hic_array=()
		qc_array=()
		qc_filename=$feather_output/$dataset_name".feather.qc"
		for i in `seq $number_of_datasets`
		do
			ds_name="dataset$i"
			hic_name="$ds_name/$hic_dir"
			#printf "$dataset1"
			#printf "$ds_name\n"
			eval ds=\$$ds_name
			eval hic=\$$hic_name
			qc_file=$(ls $ds/*.qc.tsv)
			awk 'FNR == 3 {rmdup=$2}; FNR == 4 {intra=$2}; {a[FNR]=$1; b[FNR]=$2; c[FNR]=$3} END{for (i=1; i<=FNR; i++) if (i != 12 && i != 13) {print a[i], b[i], c[i]} else {if( i== 12) {print a[i], b[i]*rmdup, c[i]}else{print a[i], b[i]*intra, c[i]}}}' $qc_file > $qc_filename".s"$i
			#echo $qc_file >> $qc_filename".s"$i
			#printf "$ds\n"
			ds_array+=($ds)
			hic_array+=($hic)
			qc_array+=($qc_filename".s"$i)
			#printf "$i\n"
			#printf '%s\n' "${ds_array[@]}"
			#printf '%s\n' "${hic_array[@]}"
			#printf '%s\n' "${qc_array[@]}"
		done
		#printf '%s\n' "${qc_array[@]}"
		$cwd/feather/concat_bedfiles.sh $feather_output $dataset_name "${ds_array[@]}"
		awk '{a[FNR]=$1; b[FNR]+=$2; c[FNR]=$3} END{for (i=1; i<=FNR; i++) print a[i], b[i], c[i]}' "${qc_array[@]}" > $qc_filename"_tmp"
		awk 'FNR == 3 {rmdup=$2}; FNR == 4 {intra=$2}; {a[FNR]=$1; b[FNR]=$2; c[FNR]=$3} END{for (i=1; i<=FNR; i++) if (i != 12 && i != 13) {print a[i], b[i], c[i]} else {if( i == 12) {print a[i], b[i]/rmdup, c[i]}else{print a[i], b[i]/intra, c[i]}}}' $qc_filename"_tmp" > $qc_filename".tsv"
		sed -i 's/ /\t/g' $qc_filename".tsv"
		if [ $generate_hic -eq 1 ]; then
			echo "$feather_output/$hic_dir"
			mkdir -p $feather_output"/"$hic_dir
			$cwd/feather/concat_hic.sh  $feather_output $dataset_name $hic_dir "${hic_array[@]}"
		fi
	else
		$python_path $cwd/feather/feather_pipe preprocess -o $feather_output -p $dataset_name -f1 $fastq1 -f2 $fastq2 -b $bwa_index -q $mapq -l $length_cutoff -t $threads -c $per_chr -j $generate_hic -a $macs2_filepath
		qc_filename=$feather_output/$dataset_name".feather.qc"
		temp_qc_file=$feather_output/tempfiles/$dataset_name".feather.qc.modified"
		#printf "dataset name:\t"$dataset_name"\n" >> $qc_filename
		#printf "MACS2 file:\t"$macs2_filename"\n" >> $qc_filename
		sed -r 's/  +/\t/g' $qc_filename > $temp_qc_file
		sed -r -i 's/ /\_/g' $temp_qc_file
		cut -f 1-2 $temp_qc_file > $temp_qc_file".cut"
		sed -i 's/\_$//g' $temp_qc_file".cut"
		paste -d"\t" $cwd/feather/qc_template.txt $temp_qc_file".cut" > $temp_qc_file".cut.tmp"
		awk '{print $1,$3,$2}' $temp_qc_file".cut.tmp" >> $qc_filename".tsv"
		sed -i 's/ /\t/g' $qc_filename".tsv"
		#awk '{printf("%s\t", $1)}' $temp_qc_file".cut" > $qc_filename".tsv"
		#printf "\n" >> $qc_filename".tsv"
		#awk '{printf("%s\t", $2)}' $temp_qc_file".cut" >> $qc_filename".tsv"
	fi
	cp "$(readlink -f $0)" $feather_output"/execution_script_copy"
	chmod 777 $feather_output
	ln -sfn $feather_output $feather_output_symlink
else
	feather_output=$(cd -P "$feather_output_symlink" && pwd)
fi

if [ $maps -eq 1 ]; then
	mkdir -p $maps_output
	echo "$dataset_name $maps_output $macs2_filepath $genomic_feat_filepath $long_bedpe_dir $short_bed_dir $bin_size $chr_count $maps_output"
	$python_path $cwd/MAPS/make_maps_runfile.py $dataset_name $maps_output $macs2_filepath $genomic_feat_filepath $long_bedpe_dir $short_bed_dir $bin_size $chr_count $maps_output $sex_chroms_to_process --BINNING_RANGE $binning_range
	echo "first"
	$python_path $cwd/MAPS/MAPS.py $maps_output"/maps_"$dataset_name".maps"
	echo "second"
	$Rscript_path $cwd/MAPS/MAPS_regression_and_peak_caller.r $maps_output $dataset_name"."$resolution"k" $bin_size $chr_count$sex_chroms $filter_file $model
	$Rscript_path $cwd/MAPS/MAPS_peak_formatting.r $maps_output $dataset_name"."$resolution"k" $fdr $bin_size
	echo "third"
	cp "$(readlink -f $0)" $maps_output"/execution_script_copy"
	chmod 777 $maps_output
	ln -sfn $maps_output $maps_output_symlink
else
	maps_output=$(cd -P "$maps_output_symlink" && pwd)
fi

############################################################################
###               Arima Genomics Post-processing and QC                  ###
############################################################################
if [ $Arima_QC -ne 1 ]; then
	echo "MAPS pipeline finished!"
	exit 0
fi

echo -e "\nPost-processing by Arima Genomics ..."

loop_file=$maps_output"/"$dataset_name"."$resolution"k."$fdr".peaks.bedpe"

if [ -f "$loop_file" ]; then
	loop_bases_file_bed=$maps_output"/"$dataset_name"."$resolution"k."$fdr".loopBases.bed"
	loop_bases_file_bed_uniq=$maps_output"/"$dataset_name"."$resolution"k."$fdr".uniqLoopBases.bed"
	cut -f1,2,3 $loop_file | grep -v "start1" > $loop_bases_file_bed
	cut -f4,5,6 $loop_file | grep -v "start2" >> $loop_bases_file_bed
	cat $loop_bases_file_bed | sort -k1,1 -k2,2n -k3,3n | uniq > $loop_bases_file_bed_uniq
	total_loops=`grep -v "start1" $loop_file | wc -l`
	#peaks_in_loops=`bedtools intersect -wa -u -a $macs2_filepath -b $loop_bases_file_bed_uniq | wc -l`
else
	total_loops=0
fi

total_peaks=`wc -l $macs2_filepath | awk '{print $1}'`
#peaks_in_loops_p=`echo "scale=4; 100 * $peaks_in_loops / $total_peaks" | bc | awk '{ printf("%.1f", $0) }'`

if [ $plot -eq 1 ]; then

	plots_dir=$outdir"/arcplot_and_metaplot/"
	mkdir $plots_dir
	echo "Generating arcplot files for WashU Epigenome Browser ..."
	if [ ! -f $loop_file ]; then
		echo -e "No bedpe loop file found, skip generating arcplot files! \n"
	else
		Rscript $cwd/bedpe2tabix.R -i $loop_file -t $plots_dir"/"$dataset_name"."$resolution"k."$fdr".arcplot"
		echo -e "Finished making arcplot files! \n"
	fi

	echo "Generating metaplot ..."
	bam_file=$feather_output"/"$dataset_name".paired.srtn.rmdup.bam"
	bam_file_sorted=$feather_output"/"$dataset_name".paired.rmdup.bam"
	bigwig_file=$plots_dir"/"$dataset_name".coverage.bigwig"
	matrix_file=$plots_dir"/"$dataset_name".coverage_matrix.tab.gz"
	metaplot=$plots_dir"/"$dataset_name".metaplot.pdf"
	heatmap=$plots_dir"/"$dataset_name".heatmap.pdf"

	echo "Creating index for the bam file ..."
	samtools index $bam_file_sorted

	echo "Generating coverage bigwig file from sequencing reads ..."
	bamCoverage \
	--bam $bam_file_sorted \
	--outFileName $bigwig_file \
	--outFileFormat bigwig \
	--ignoreDuplicates \
	--numberOfProcessors max

	# Generate a matrix file with +/- 10kb from the center of each region in bed file
	echo "Generating coverage matrix file ..."
	computeMatrix reference-point \
	--scoreFileName $bigwig_file \
	--regionsFileName $macs2_filepath \
	--referencePoint center \
	--beforeRegionStartLength 10000 \
	--afterRegionStartLength 10000 \
	--outFileName $matrix_file \
	--numberOfProcessors max

	# Use plotProfile to plot a line graph of the average of the matrix
	echo "Generating metaplot image ..."
	plotProfile \
	--matrixFile $matrix_file \
	--outFileName $metaplot \
	--yAxisLabel "Coverage" \
	--refPointLabel "peak center" \
	--regionsLabel $dataset_name \
	--samplesLabel "HiChIP Signal Enrichment at ChIP Peaks" \
	--yMin 0

	# Use plotHeatmap to plot a heatmap of the average of the matrix
	echo "Generating heatmap image ..."
	plotHeatmap \
	--matrixFile $matrix_file \
	--outFileName $heatmap \
	--xAxisLabel "" \
	--yAxisLabel "Coverage" \
	--refPointLabel "peak center" \
	--regionsLabel $dataset_name \
	--samplesLabel "HiChIP Signal Enrichment at ChIP Peaks" \
	--colorList "white,darkblue" \
	--heatmapHeight 12 \
	--yMin 0

	echo -e "Finished making metaplot and heatmap iamges! \n"

fi

feather_qc=$feather_output"/"$dataset_name".feather.qc"
maps_qc=$maps_output"/"$dataset_name".maps.qc"

if [ ! -f "$feather_qc" ] || [ ! -f "$maps_qc" ]; then
  	echo "*.qc file does not exist! Exiting ..."
	exit 1
fi

#echo peaks_in_loops=$peaks_in_loops
#echo peaks_in_loops_p=$peaks_in_loops_p

# From MAPS, please note the MAPS version!
pairs_raw=`grep "number of sequencing pairs" $feather_qc | awk '{print $5}'`
#pairs_uniq_mapped_MAPQ_ge_30=`grep "number of uniquely mapped pairs (MAPQ >= 30)" $feather_qc | awk '{print $9}'`
pairs_uniq_mapped_MAPQ_ge_30=`grep "number of un" $feather_qc | awk '{print $9}'`
pairs_uniq_mapped_MAPQ_ge_30_p=`echo "scale=4; 100 * $pairs_uniq_mapped_MAPQ_ge_30 / $pairs_raw" | bc | awk '{ printf("%.1f", $0) }'`
pairs_deduped=`grep "number of pairs after duplicate removal" $feather_qc | awk '{print $7}'`
pairs_deduped_p=`echo "scale=4; 100 * $pairs_deduped / $pairs_raw" | bc | awk '{ printf("%.1f", $0) }'`

#echo pairs_raw=$pairs_raw
#echo pairs_uniq_mapped_MAPQ_ge_30=$pairs_uniq_mapped_MAPQ_ge_30
#echo pairs_uniq_mapped_MAPQ_ge_30_p=$pairs_uniq_mapped_MAPQ_ge_30_p
#echo pairs_deduped=$pairs_deduped
#echo pairs_deduped_p=$pairs_deduped_p

# From MAPS, please note the MAPS version!
pairs_intra=`grep "number of intra" $feather_qc | awk '{print $5}'`
reads_intra=$(($pairs_intra * 2))
pairs_inter=`grep "number of inter" $feather_qc | awk '{print $5}'`
reads_inter=$(($pairs_inter * 2))

#echo pairs_intra=$pairs_intra
#echo pairs_inter=$pairs_inter

# HiChIP Data Features
short_VIPs=`wc -l $feather_output"/"$dataset_name".shrt.vip.bed" | awk '{print $1}'`
short_VIPs_file=$feather_output"/"$dataset_name".shrt.vip.bed"
short_VIPs_file_sorted=$feather_output"/"$dataset_name".shrt.vip.sort.bed"
cut -f1-3 $short_VIPs_file | sort -k1,1 -k2,2n -k3,3n > $short_VIPs_file_sorted
short_VIPs_with_peaks=`bedtools intersect -a $short_VIPs_file_sorted -b $macs2_filepath | wc -l`
short_VIPs_with_peaks_p=`echo "scale=4; 100 * $short_VIPs_with_peaks / $short_VIPs" | bc | awk '{ printf("%.1f", $0) }'`

AND=`grep AND $maps_qc | awk '{sum=sum+$2}; END {print sum}'`
XOR=`grep XOR $maps_qc | awk '{sum=sum+$2}; END {print sum}'`
NOT=`grep NOT $maps_qc | awk '{sum=sum+$2}; END {print sum}'`

########### New FRIPs calculation ############
long_intra_file=$feather_output"/"$dataset_name".long.intra.bedpe"
long_intra_filtered_file=$feather_output"/"$dataset_name".long.intra.filtered.bedpe"

cat $long_intra_file | awk -v resolution=$resolution 'BEGIN {resolution=resolution * 1000 } {R1_midpoint=(($2+$3)/2); R1_bin=int(R1_midpoint / resolution); R1_bin_start=(R1_bin * resolution); R1_bin_end=((R1_bin * resolution) + resolution - 1); R1_bin_mid=((R1_bin_start + R1_bin_end)/2); R2_midpoint=(($5+$6)/2); R2_bin=int(R2_midpoint / resolution); R2_bin_start=(R2_bin * resolution); R2_bin_end=((R2_bin * resolution) + resolution - 1); R2_bin_mid=((R2_bin_start + R2_bin_end)/2); distance=(sqrt((R1_bin_mid - R2_bin_mid)^2)); {if(distance >= 10000 && distance <= 2000000) { print $1"\t"R1_bin_start"\t"R1_bin_end"\t"$4"\t"R2_bin_start"\t"R2_bin_end } } }' | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n > $long_intra_filtered_file

# Calculate Arima's AND, XOR and NOT
ALL_Arima=`cat $long_intra_filtered_file | wc -l`
AND_Arima=`cut -f4,15 $maps_output"/reg_raw.chr"*"."*"."$resolution"k.and" | grep -v "dist" | awk '{ if($2 > 1) a+=$1 } END {print a}'`
XOR_Arima=`cut -f4,15 $maps_output"/reg_raw.chr"*"."*"."$resolution"k.xor" | grep -v "dist" | awk '{ if($2 > 1) a+=$1 } END {print a}'`
FRIPs=$(( $AND_Arima + $XOR_Arima ))
NOT_Arima=$(( $ALL_Arima - $FRIPs ))
FRIPs_p=`echo "scale=4; 100 * $FRIPs / $ALL_Arima" | bc | awk '{ printf("%.1f", $0) }'`

bins_with_peaks=`awk '{ bin1 = int($2 / 5000); bin2 = int($3 / 5000); for(i = bin1; i <= bin2; i++) print $1"_"i }' $macs2_filepath | sort | uniq | wc -l`
FRIPs_per_peak_bins=`echo "scale=4; $FRIPs / $bins_with_peaks" | bc | awk '{ printf("%.1f", $0) }'`
FRIPs_per_peak_bins_per_million_deduped=`echo "scale=5; ($FRIPs / $bins_with_peaks) * 1000000 / $pairs_deduped" | bc | awk '{ printf("%.2f", $0) }'`
target_FRIPs_per_peak_bin=382
mapping_rate="0.70"
sequence_uniqueness_rate="0.60"
rounding_level="10000000"
target_raw_pairs=`echo "scale=6; ( 1000000 * $target_FRIPs_per_peak_bin / $FRIPs_per_peak_bins_per_million_deduped / ($mapping_rate * $sequence_uniqueness_rate) )" | bc | awk -v rounding_level=$rounding_level '{rounded_num = ((int($0 / rounding_level) + 1) * rounding_level)} { printf("%d", rounded_num) }'`

#echo short_VIPs=$short_VIPs
#echo short_VIPs_with_peaks=$short_VIPs_with_peaks
#echo short_VIPs_with_peaks_p=$short_VIPs_with_peaks_p
#echo AND=$AND
#echo XOR=$XOR
#echo NOT=$NOT
#echo ALL_Arima=$ALL_Arima
#echo AND_Arima=$AND_Arima
#echo XOR_Arima=$XOR_Arima
#echo NOT_Arima=$NOT_Arima
#echo FRIPs=$FRIPs
#echo FRIPs_p=$FRIPs_p
#echo bins_with_peaks=$bins_with_peaks
#echo FRIPs_per_peak_bins=$FRIPs_per_peak_bins
#echo FRIPs_per_peak_bins_per_million_deduped=$FRIPs_per_peak_bins_per_million_deduped
#echo target_raw_pairs=$target_raw_pairs

echo -e "Removing duplicates using Picard ...\n"
MAX_memory=$(($threads * 4 - 8))
INPUT=$feather_output"/"$dataset_name".paired.srt.bam" #JMB: more descriptive variable name
OUTPUT=$feather_output"/"$dataset_name".paired.nodup.bam" #JMB: more descriptive variable name
METRICS_FILE=$feather_output"/"$dataset_name".paired.nodup.metrics.txt" #JMB: more descriptive variable name
$java_path -Xmx${MAX_memory}G -XX:-UseGCOverheadLimit -Djava.io.tmpdir=temp/ -jar $picard_path MarkDuplicates INPUT=$INPUT OUTPUT=$OUTPUT METRICS_FILE=$METRICS_FILE TMP_DIR=temp/ ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=TRUE

if [ ! -f $METRICS_FILE ]; then
	echo "Could not find metrics file generated by MarkDuplicates! Exiting ..."
	exit 1
else
	echo -e "\nDone removing duplicates! \n"
fi

# From picard.jar MarkDuplicates
pairs_mapped=$( awk 'FNR==8 {print}' $METRICS_FILE | cut -f3 )
pairs_duplicates=$( awk 'FNR==8 {print}' $METRICS_FILE | cut -f7 )
pairs_opt_duplicates=$( awk 'FNR==8 {print}' $METRICS_FILE | cut -f8 )
pairs_PCR_duplicates=$(($pairs_duplicates - $pairs_opt_duplicates))
pairs_PCR_duplicates_p=`echo "scale=4; 100 * $pairs_PCR_duplicates / $pairs_mapped" | bc | awk '{ printf("%.1f", $0) }'`

#echo pairs_mapped=$pairs_mapped
#echo pairs_duplicates=$pairs_duplicates
#echo pairs_PCR_duplicates=$pairs_PCR_duplicates
#echo pairs_PCR_duplicates_p=$pairs_PCR_duplicates_p

tmp_file=$feather_output"/cov.tmp"
reads_all=`samtools view $OUTPUT | wc -l`
pairs_all=$(($reads_all / 2))
samtools view $OUTPUT | awk '{ if($7=="=") print $9 }' > $tmp_file
reads_intra=`wc -l $tmp_file | awk '{print $1}'`
pairs_intra=$(($reads_intra / 2))
pairs_intra_p=`echo "scale=4; 100 * $pairs_intra / $pairs_all" | bc | awk '{ printf("%.1f", $0) }'`
reads_inter=$(($reads_all - $reads_intra))
pairs_inter=$(($reads_inter / 2))
pairs_inter_p=`echo "scale=4; 100 * $pairs_inter / $pairs_all" | bc | awk '{ printf("%.1f", $0) }'`
pairs_intra_ge_15kb=`awk '{if($1 >= 15000) a++} END {print a}' $tmp_file`
pairs_intra_ge_15kb_p=`echo "scale=4; 100 * $pairs_intra_ge_15kb / $pairs_all" | bc | awk '{ printf("%.1f", $0) }'`

#echo pairs_all=$pairs_all
#echo pairs_intra=$pairs_intra
#echo pairs_intra_p=$pairs_intra_p
#echo pairs_inter=$pairs_inter
#echo pairs_inter_p=$pairs_inter_p

QC_result_deep=$outdir"/"Arima_QC_deep.txt
QC_result_shallow=$outdir"/"Arima_QC_shallow.txt

#JMB: should we print out "AND", "XOR", and "NOT" or should we print "AND_Arima", "XOR_Arima" and "NOT_Arima"?

header_deep=("Sample Name" "Raw PE Reads" "Mapped Reads" "% Mapped Reads" "PCR Dups" "% PCR Dups" "Mapped and De-duped" "% Mapped and De-duped" "Loops" "LR FRIPS per Peak Bin per Million Deduped PE Reads" "LR FRIPS per Peak Bin" "LR FRIPs" "% LR FRIPs" "AND pairs" "XOR pairs" "NOT pairs" "INTRA pairs" "% INTRA pairs" "INTRA > 15kb pairs" "% INTRA > 15kb pairs" "INTER pairs" "% INTER pairs" "ChIP Peaks" "Peak Bins" "Short VIPs" "Short VIPs in Peaks" "% Short VIPs in Peaks")
IFS=$'\t'; echo "${header_deep[*]}" > $QC_result_deep

result_deep=($dataset_name $pairs_raw $pairs_uniq_mapped_MAPQ_ge_30 $pairs_uniq_mapped_MAPQ_ge_30_p $pairs_PCR_duplicates $pairs_PCR_duplicates_p $pairs_deduped $pairs_deduped_p $total_loops $FRIPs_per_peak_bins_per_million_deduped $FRIPs_per_peak_bins $FRIPs $FRIPs_p $AND_Arima $XOR_Arima $NOT_Arima $pairs_intra $pairs_intra_p $pairs_intra_ge_15kb $pairs_intra_ge_15kb_p $pairs_inter $pairs_inter_p $total_peaks $bins_with_peaks $short_VIPs $short_VIPs_with_peaks $short_VIPs_with_peaks_p)
IFS=$'\t'; echo "${result_deep[*]}" >> $QC_result_deep

header_shallow=("Sample Name" "Raw PE Reads" "Mapped Reads" "% Mapped Reads" "PCR Dups" "% PCR Dups" "Mapped and De-duped" "% Mapped and De-duped" "Target Raw PE Reads" "LR FRIPS per Peak Bin per Million Deduped PE Reads" "LR FRIPS per Peak Bin" "LR FRIPs" "% LR FRIPs" "AND pairs" "XOR pairs" "NOT pairs" "INTRA pairs" "% INTRA pairs" "INTRA > 15kb pairs" "% INTRA > 15kb pairs" "INTER pairs" "% INTER pairs" "ChIP Peaks" "Peak Bins" "Short VIPs" "Short VIPs in Peaks" "% Short VIPs in Peaks")
IFS=$'\t'; echo "${header_shallow[*]}" > $QC_result_shallow

result_shallow=($dataset_name $pairs_raw $pairs_uniq_mapped_MAPQ_ge_30 $pairs_uniq_mapped_MAPQ_ge_30_p $pairs_PCR_duplicates $pairs_PCR_duplicates_p $pairs_deduped $pairs_deduped_p $target_raw_pairs $FRIPs_per_peak_bins_per_million_deduped $FRIPs_per_peak_bins $FRIPs $FRIPs_p $AND_Arima $XOR_Arima $NOT_Arima $pairs_intra $pairs_intra_p $pairs_intra_ge_15kb $pairs_intra_ge_15kb_p $pairs_inter $pairs_inter_p $total_peaks $bins_with_peaks $short_VIPs $short_VIPs_with_peaks $short_VIPs_with_peaks_p)
IFS=$'\t'; echo "${result_shallow[*]}" >> $QC_result_shallow

echo "MAPS and Arima QC pipeline finished successfully!"
echo "Please download the QC result from: $QC_result_deep $QC_result_shallow and then copy the contents to the corresponding tables in the QC worksheet."
