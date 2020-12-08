#!/bin/bash
############################################################################
###                       Arima MAPS v2 QC Metrics                       ###
############################################################################

# Use this config file to run Arima HiChIP data, created using the Arima Hi-C+ kit, through the MAPS pipeline
# This pipeline was modified from MAPS: https://github.com/ijuric/MAPS
# All modifications to this code are in post processing steps to generate Arima specific quality control metrics and do not affect the MAPS algorithm itself.
# For technical assistance please contact Ming Hu at hum@ccf.org
# Please use the following citation for referencing MAPS: Juric I, Yu M, Abnousi A, Raviram R, Fang R, Zhao Y, et al. (2019) MAPS: Model-based analysis of long-range chromatin interactions from PLAC-seq and HiChIP experiments. PLoS Comput Biol 15(4): e1006982. https://doi.org/10.1371/ journal.pcbi.1006982

############################################################################
###                       Computing Environment                          ###
############################################################################

# Please install the following dependencies and add them to your PATH variable to ensure MAPS will execute in your computing environment:

# Computing Resources: For shallow sequencing (0.5 - 2 million raw paired-end reads), the Arima MAPS pipeline requires 12 CPU cores with 48 GB RAM.  The shallow sequencing analysis should complete in less than 2 hrs.  For deep sequencing (200 – 600 million raw paired-end reads), we recommend 16 - 20 CPU cores with at least 64 - 80 GB RAM.  Samples with 200 million raw paired-end reads will run in about 48 hours with the recommended computational resources.  Additional resources can be added to decrease the analysis time.

# Dependencies:

# Python 3.4 (or later): pandas (v0.20.3 or later), numpy (v1.13.1 or later), itertools (v3.2), pysam (v0.15.2 or later), pybedtools (v0.8.0 or later), deeptools (v3.3.0 or later)
# R 3.4.3 packages: argparse (v2.0.1), VGAM (v1.1-2), data.table (v1.12.8)
# other: bedtools (v2.27.1), htslib (v1.10), samtools (v1.10), bcftools (v1.10).

# see https://github.com/ijuric/MAPS/tree/master/Arima_Genomics/README.md for installation help

# please provide the path to the following tools:

#python_path=/home/xiangz/tools/anaconda3/bin/python
# Example: python_path=~/anaconda3/bin/python
#Rscript_path=/opt/R/bin/Rscript
# Example: Rscript_path=/opt/R/bin/Rscript
cwd=$(dirname $0)
# cwd=/home/xiangz/tools/MAPS/bin/
# Example: cwd=~/tools/MAPS/bin/
#MACS2_path=/home/xiangz/tools/anaconda3/bin/macs2
# Can be set to empty ("") if you provide your own ChIP peak file

############################################################################
###                          Customer Inputs                             ###
############################################################################

call_peaks=0 # "1" to call peaks using MACS2, "0" to use ChIP peak file provided
peak_type="" # Required if call_peaks=1. Type of peaks to call. Must be either "broad" (H3K4me3) or "narrow" (CTCF)
feather=1 # "1" to run feather, "0" to skip
maps=1 # "1" to run MAPS on data previously processed with feather, "0" to skip
# setting both feather and MAPS to zero will run the Arima QC metrics only on previous processed data
number_of_datasets=1 # number of biological replicas. Can be 1,2, … N. If set to 1, script assumes that user is running MAPS for a particular biological replica. If number_of_datasets is larger than 1, script assumes that user is running MAPS on merged data, and it will look for .bedpe and .bed files of each replica defined in dataset1, dataset2,..., datasetN lines (See Additional Parameters section)

# Removed in v2.0
# dataset_name="Arima_MAPS_test" # sample name with no file extention, Ex: "sample1"
# fastq_format=".fastq.gz" # format that MAPS should expect the data in. Ex: "_001.fastq.gz"
# pair_format="_" # fastq file names must be named "08-H025-08_S13_L001"[pair_format][R1 or R2][fastq_format]. Ex: "_" or "."
# fastq_dir="/oasis/tscc/scratch/xiangz/MAPS/v2/Arima_FTP/fastq/" # absolute path to the directory the data is in. Ex: "/home/sample1/fastq/"

# Added in v2.0
fastq_dir_with_file_prefix="" # absolute path and file prefix of the fastq files, up to "_R1" or ".R1"

outdir="" # directory which the feather and MAPS outputs will be placed
macs2_filepath="" # path to the ChIP peak file that will be used for loop calling with MAPS. Required if call_peaks=0. Leave it empty ("") if call_peaks=1
organism="" # organism of the genomic feature file to be used, options: "mm9", "mm10", "hg19", and "hg38"
bwa_index="" # absolute path to the reference genome sequence (.fa) which will be used to derive the BWA index as well. Ex: "/home/reference_sequence/hg19.fa"
threads=8 # number of threads for MAPS to run on. Use 4-8 for shallow sequencing and 12-20 for deep sequencing.
patterned_flowcell="" # Use "1" for deep sequencing and "0" for shallow sequencing datasets. "1" if the data was sequenced on a patterned flowcell, "0" if non-patterned flowcell. This is used for calculating optical duplicates and PCR duplicate rates. # v1.9_update

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
# number_of_datasets > 1 (merging existing Feather processed datasets)
# specify as many datasets as required

#dataset1=""
#dataset2=""
#dataset3=""
#dataset4=""
#...

###SET THESE VARIABLES ONLY IF FEATHER = 0 AND YOU WANT TO RUN
###USING A SPECIFIC FEATHER OUTPUT RATHER THAN $datasetname_Current

feather_output_symlink=""


############################################################################
###                      Command Line Arguments                          ###
############################################################################

usageHelp="Usage: ${0##*/} [-C call_peaks] [-p peak_type] [-F feather] [-M maps]
       [-I fastq_dir_with_file_prefix] [-O outdir] [-m macs2_filepath]
       [-o organism] [-b bwa_index] [-t threads] [-f patterned_flowcell]
       [-P plot] [-H generate_hic] [-s bin_size] [-r binning_range] [-d fdr]
       [-Q mapq] [-l length_cutoff] [-h] \n"
call_peaksHelp="* [-C call_peaks]: 0 (default) to use ChIP peak file provided, 1 to call peaks
    using MACS2"
peak_typeHelp="* [-p peak_type]: broadness of chromatin factor. Required if call_peaks=1. Must
    be either \"broad\" (H3K4me3) or \"narrow\" (CTCF). Both choices result in broad
    peaks being called by MACS2."
featherHelp="* [-F feather]: 1 (default) to run feather, 0 to skip"
mapsHelp="* [-M maps]: 1 (default) to run MAPS on data processed with feather, 0 to skip"
fastq_dir_with_file_prefixHelp="* [-I fastq_dir_with_file_prefix]: absolute path and file prefix of the fastq
    files, up to \"_R1\" or \".R1\""
outdirHelp="* [-O outdir]: directory which the feather and MAPS outputs will be placed"
macs2_filepathHelp="* [-m macs2_filepath]: path to the ChIP peak file that will be used for loop
    calling with MAPS. Required if call_peaks=0. Ignore or leave it empty (\"\")
    if call_peaks=1"
organismHelp="* [-o organism]: organism of the genomic feature file to be used, options:
    \"mm9\", \"mm10\", \"hg19\", and \"hg38\""
bwa_indexHelp="* [-b bwa_index]: absolute path to the reference genome sequence (.fa) which will
    be used to derive the BWA index as well. Ex: \"/home/reference_sequence/hg19.fa\""
threadsHelp="* [-t threads]: number of threads for MAPS to run on. Use 4-8 for shallow sequencing
    and 12-20 for deep sequencing."
patterned_flowcellHelp="* [-f patterned_flowcell]: Use 1 for deep sequencing and 0 for shallow sequencing
    datasets. \"1\" if the data was sequenced on a patterned flowcell, \"0\" if
    non-patterned flowcell. This is used for calculating optical duplicates and
    PCR duplicate rates."
plotHelp="* [-P plot]: 1 (default) to generate plots, 0 to skip"
generate_hicHelp="* [-H generate_hic]: 1 (default) to generate a .hic file for visualization with
    Juicer, 0 to skip"
bin_sizeHelp="* [-s bin_size]: resolution of the loops called. Default: 5000"
binning_rangeHelp="* [-r binning_range]: Maximum distance for loop calling. Default: 2000000"
fdrHelp="* [-d fdr]: Do not change! False Discovery Rate threshold: 1 = 0.1,
    2 (default) = 0.01, 3 = 0.001, ect..."
mapqHelp="* [-Q mapq]: Phred scaled mapping quality threshold. Default: 30"
length_cutoffHelp="* [-l length_cutoff]: Minimum genomic distance. Default: 1000"
helpHelp="* [-h]: print this help and exit"

printHelpAndExit() {
    echo -e "$usageHelp"
	echo -e "$call_peaksHelp"
    echo -e "$peak_typeHelp"
    echo -e "$featherHelp"
    echo -e "$mapsHelp"
    echo -e "$fastq_dir_with_file_prefixHelp"
    echo -e "$outdirHelp"
    echo -e "$macs2_filepathHelp"
    echo -e "$organismHelp"
    echo -e "$bwa_indexHelp"
    echo -e "$threadsHelp"
    echo -e "$patterned_flowcellHelp"
    echo -e "$plotHelp"
    echo -e "$generate_hicHelp"
    echo -e "$bin_sizeHelp"
    echo -e "$binning_rangeHelp"
    echo -e "$fdrHelp"
    echo -e "$mapqHelp"
    echo -e "$length_cutoffHelp"
    echo "$helpHelp"
    exit "$1"
}

while getopts "b:C:d:f:F:hH:I:l:m:M:o:O:p:P:Q:r:s:t:" opt; do
    case $opt in
	h) printHelpAndExit 0;;
	C) call_peaks=$OPTARG ;;
	p) peak_type=$OPTARG ;;
	F) feather=$OPTARG ;;
	M) maps=$OPTARG ;;
	I) fastq_dir_with_file_prefix=$OPTARG ;;
	O) outdir=$OPTARG ;;
	m) macs2_filepath=$OPTARG ;;
	o) organism=$OPTARG ;;
	b) bwa_index=$OPTARG ;;
	t) threads=$OPTARG ;;
	f) patterned_flowcell=$OPTARG ;;
	P) plot=$OPTARG ;;
	H) generate_hic=$OPTARG ;;
	s) bin_size=$OPTARG ;;
	r) binning_range=$OPTARG ;;
	d) fdr=$OPTARG ;;
	Q) mapq=$OPTARG ;;
	l) length_cutoff=$OPTARG ;;
	[?]) printHelpAndExit 1;;
    esac
done

# Sanity checks
if [ -z "$fastq_dir_with_file_prefix" ]; then
    echo "Please provide input files (-I)!"
    printHelpAndExit 1
fi

if [ -z "$outdir" ]; then
    echo "Please provide an output folder (-O)!"
    printHelpAndExit 1
fi

if [ -z "$organism" ]; then
    echo "Please provide an organism (-o)!"
    printHelpAndExit 1
fi

if [ ! -f "$bwa_index" ]; then
    echo "Please provide a reference genome file (-b)!"
    printHelpAndExit 1
fi

if [ -z "$patterned_flowcell" ]; then
    echo "Please provide a value for patterned flowcell (-f)!"
    printHelpAndExit 1
fi

echo "User Defined Inputs:"
echo call_peaks=$call_peaks
if [ "$call_peaks" -eq 1 ]; then
    if [ -z "$peak_type" ]; then
        echo "Broadness of chromatin factor (-p) is required when call_peaks=1!"
        printHelpAndExit 1
    else
        echo peak_type=$peak_type
    fi
fi
echo feather=$feather
echo maps=$maps
echo fastq_dir_with_file_prefix=$fastq_dir_with_file_prefix
echo outdir=$outdir
if [ "$call_peaks" -eq 0 ]; then
	if [ ! -f "$macs2_filepath" ]; then
		echo "ChIP peak file (-m) is required for loop calling when call_peaks=0!"
		printHelpAndExit 1
	else
		echo macs2_filepath=$macs2_filepath
	fi
fi
echo organism=$organism
echo bwa_index=$bwa_index
echo threads=$threads
echo patterned_flowcell=$patterned_flowcell
echo plot=$plot
echo generate_hic=$generate_hic
echo bin_size=$bin_size
echo binning_range=$binning_range
echo fdr=$fdr
echo mapq=$mapq
echo length_cutoff=$length_cutoff
echo

# You can specify the absolute path for genomics feature files below
MAPS_DIR=`dirname $cwd`
if [ $organism == "mm10" ]; then
	genomic_feat_filepath=$MAPS_DIR"/Arima_Genomics/genomic_features/mm10_F_GC_M_arima_5Kb_el.txt"
	chr_count=19
elif [ $organism == "mm9" ]; then
	genomic_feat_filepath=$MAPS_DIR"/Arima_Genomics/genomic_features/mm9_F_GC_M_arima_5Kb_el.txt"
	chr_count=19
elif [ $organism == "hg19" ]; then
	genomic_feat_filepath=$MAPS_DIR"/Arima_Genomics/genomic_features/hg19_F_GC_M_arima_5Kb_el.txt"
	chr_count=22
elif [ $organism == "hg38" ]; then
	genomic_feat_filepath=$MAPS_DIR"/Arima_Genomics/genomic_features/hg38_F_GC_M_arima_5Kb_el.txt"
	chr_count=22
fi

optical_duplicate_distance=$([ $patterned_flowcell == 1 ] && echo 2500 || echo 100) # v1.9_update


############################################################################
###                           MAPS Pipeline                              ###
############################################################################

DATE=`date '+%Y%m%d_%H%M%S'`

# Removed in v2.0
#fastq1=$fastq_dir/$dataset_name$pair_format"R1"$fastq_format
#fastq2=$fastq_dir/$dataset_name$pair_format"R2"$fastq_format

# Added in v2.0
fastq1=$fastq_dir_with_file_prefix"[._]R1*"
fastq2=$fastq_dir_with_file_prefix"[._]R2*"
fastq_dir=$(dirname $fastq_dir_with_file_prefix)
dataset_name=$(basename $fastq_dir_with_file_prefix)
fastq_format=$(echo $fastq1 | sed 's/^.*[._]R1\(.*$\)/\1/')

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

### Added support for multiple samples
reads_all_duplicates=0
reads_opt_duplicates=0

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

### Added support for multiple samples
			bam_rmdup=$(ls $ds/*.paired.rmdup.bam)
			bam_array+=($bam_rmdup)
			dup_stats_file=$(ls $ds/*.paired.fixmated.markdup.stats)

			reads_all_dups[$i]=$(cat $dup_stats_file | grep -P "DUPLICATE PAIR:" | cut -d' ' -f3)
			reads_opt_dups[$i]=$(cat $dup_stats_file | grep -P "DUPLICATE PAIR OPTICAL:" | cut -d' ' -f4)

			reads_all_duplicates=$(( reads_all_duplicates + reads_all_dups[$i] ))
			reads_opt_duplicates=$(( reads_opt_duplicates + reads_opt_dups[$i] ))

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

### Added support for multiple samples
		merged_bam=$feather_output"/"$dataset_name".paired.rmdup.bam"
		samtools merge $merged_bam -@ $threads  "${bam_array[@]}"

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


		############################################################################
		###                           Peak Calling                               ###
		############################################################################
		if [ $call_peaks -eq 1 ]; then
			if [[ $organism == "hg19" || $organism == "hg38" ]]; then
				genome="hs"
			elif [[ $organism == "mm9" || $organism == "mm10" ]]; then
				genome="mm"
			fi

			python $cwd/feather/feather_pipe preprocess -o $feather_output -p $dataset_name -f1 $fastq1 -f2 $fastq2 -b $bwa_index -q $mapq -l $length_cutoff -t $threads -c $per_chr -j $generate_hic -d $optical_duplicate_distance

			echo -e "\nCalling $peak_type peaks using MACS2 ..."
			shortVIP_BAM=$feather_output/tempfiles/$dataset_name".shrt.bam"
			if [ $peak_type == "broad" ]; then
				macs2 callpeak -t $shortVIP_BAM -n $dataset_name -g $genome --broad --nolambda --broad-cutoff 0.3 --outdir $outdir/MACS2_peaks/
			elif [ $peak_type == "narrow" ]; then
				macs2 callpeak -t $shortVIP_BAM -n $dataset_name -g $genome --broad --broad-cutoff 0.2 --outdir $outdir/MACS2_peaks/
			fi
			rm $outdir/MACS2_peaks/${dataset_name}_peaks.gappedPeak
			macs2_filepath=$outdir/MACS2_peaks/${dataset_name}_peaks.broadPeak

			echo -e "Finished calling peaks using MACS2. The output peak file is:"
			echo -e "$macs2_filepath\n"
		else
			python $cwd/feather/feather_pipe preprocess -o $feather_output -p $dataset_name -f1 $fastq1 -f2 $fastq2 -b $bwa_index -q $mapq -l $length_cutoff -t $threads -c $per_chr -j $generate_hic -a $macs2_filepath -d $optical_duplicate_distance
		fi


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
	python $cwd/MAPS/make_maps_runfile.py $dataset_name $maps_output $macs2_filepath $genomic_feat_filepath $long_bedpe_dir $short_bed_dir $bin_size $chr_count $maps_output $sex_chroms_to_process --BINNING_RANGE $binning_range
	echo "first"
	python $cwd/MAPS/MAPS.py $maps_output"/maps_"$dataset_name".maps"
	echo "second"
	Rscript $cwd/MAPS/MAPS_regression_and_peak_caller.r $maps_output $dataset_name"."$resolution"k" $bin_size $chr_count$sex_chroms $filter_file $model
	Rscript $cwd/MAPS/MAPS_peak_formatting.r $maps_output $dataset_name"."$resolution"k" $fdr $bin_size
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

echo -e "\nPost-processing by Arima Genomics ..."

loop_file=$maps_output"/"$dataset_name"."$resolution"k."$fdr".sig3Dinteractions.bedpe"

if [ -f "$loop_file" ]; then
	total_loops=`grep -v "start1" $loop_file | wc -l`
else
	total_loops=0
fi

if [ $plot -eq 1 ]; then

	plots_dir=$outdir"/arcplot_and_metaplot/"
	[ -d $plots_dir ] || mkdir $plots_dir
	echo "Generating arcplot files for WashU Epigenome Browser ..."
	if [ ! -f $loop_file ]; then
		echo -e "No bedpe loop file found, skip generating arcplot files! \n"
	else
		Rscript $cwd/utils/bedpe2tabix.R -i $loop_file -t $plots_dir"/"$dataset_name"."$resolution"k."$fdr".arcplot"
		rm $plots_dir"/"$dataset_name"."$resolution"k."$fdr".arcplot"
		echo -e "Finished making arcplot files! \n"
	fi

	echo "Generating metaplot ..."
	bam_file_sorted=$feather_output"/"$dataset_name".paired.rmdup.bam"
	bigwig_file=$plots_dir"/"$dataset_name".coverage.bigwig"
	matrix_file=$plots_dir"/"$dataset_name".coverage_matrix.tab.gz"
	metaplot=$plots_dir"/"$dataset_name".metaplot.pdf"
	heatmap=$plots_dir"/"$dataset_name".heatmap.pdf"
	#enrichment_file=$plots_dir"/"$dataset_name".enrichment.txt"

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
	# echo "Generating metaplot image ..."
	# plotProfile \
	# --matrixFile $matrix_file \
	# --outFileName $metaplot \
	# --yAxisLabel "Coverage" \
	# --refPointLabel "peak center" \
	# --regionsLabel $dataset_name \
	# --samplesLabel "HiChIP Signal Enrichment at ChIP Peaks" \
	# --yMin 0

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

	echo "Calculating enrichment score of the metaplot ..."
	enrichment=$(zcat $matrix_file | awk -F $'\t' 'BEGIN {background = 0; peak = 0} {background = background + $7; peak = peak + $1007} END {enrichment = peak / background; printf("%.2f", enrichment)}');
	printf "%.2f\n" $enrichment
	#printf "%s\t%.2f\n" $dataset_name $enrichment > $enrichment_file

	echo -e "Finished making metaplot and heatmap images! \n"

fi

# Calculate basic data processing statistics

### Added support for multiple samples
feather_qc=$feather_output"/"$dataset_name".feather.qc.tsv"
maps_qc=$maps_output"/"$dataset_name".maps.qc"

if [ ! -f "$feather_qc" ]; then
  	echo "*.feather.qc file does not exist! Exiting ..."
	exit 1
fi

if [ ! -f "$maps_qc" ]; then
  	echo "*.maps.qc file does not exist! Exiting ..."
	exit 1
fi

### Added support for multiple samples
pairs_raw=`grep "number_of_sequencing_pairs" $feather_qc | awk '{print $2}'`
pairs_uniq_mapped_MAPQ_ge_30=`grep "number_of_uniquely_mapped_pairs_(MAPQ_>=_30)" $feather_qc | awk '{print $2}'`
pairs_uniq_mapped_MAPQ_ge_30_p=`echo "scale=4; 100 * $pairs_uniq_mapped_MAPQ_ge_30 / $pairs_raw" | bc | awk '{ printf("%.1f", $0) }'`
pairs_mapped_deduped=`grep "number_of_pairs_after_duplicate_removal" $feather_qc | awk '{print $2}'`
pairs_mapped_deduped_p=`echo "scale=4; 100 * $pairs_mapped_deduped / $pairs_raw" | bc | awk '{ printf("%.1f", $0) }'`

### Added support for multiple samples
if [ $number_of_datasets -eq 1 ]; then
	dup_stats=$feather_output"/"$dataset_name".paired.fixmated.markdup.stats"
	reads_all_duplicates=$(cat $dup_stats | grep -P "DUPLICATE PAIR:" | cut -d' ' -f3)
	reads_opt_duplicates=$(cat $dup_stats | grep -P "DUPLICATE PAIR OPTICAL:" | cut -d' ' -f4)
fi

pairs_all_duplicates=$(( $reads_all_duplicates / 2 ))
pairs_all_duplicates_p=`echo "scale=4; 100 * $pairs_all_duplicates / $pairs_uniq_mapped_MAPQ_ge_30" | bc | awk '{ printf("%.1f", $0) }'`
pairs_opt_duplicates=$(( $reads_opt_duplicates / 2 ))
pairs_opt_duplicates_p=`echo "scale=4; 100 * $pairs_opt_duplicates / $pairs_uniq_mapped_MAPQ_ge_30" | bc | awk '{ printf("%.1f", $0) }'`
pairs_PCR_duplicates=$(($pairs_all_duplicates - $pairs_opt_duplicates))
pairs_PCR_duplicates_p=`echo "scale=4; 100 * $pairs_PCR_duplicates / $pairs_uniq_mapped_MAPQ_ge_30" | bc | awk '{ printf("%.1f", $0) }'`

echo pairs_raw=$pairs_raw
echo pairs_uniq_mapped_MAPQ_ge_30=$pairs_uniq_mapped_MAPQ_ge_30
echo pairs_uniq_mapped_MAPQ_ge_30_p=$pairs_uniq_mapped_MAPQ_ge_30_p
echo pairs_mapped_deduped=$pairs_mapped_deduped
echo pairs_mapped_deduped_p=$pairs_mapped_deduped_p
echo pairs_all_duplicates=$pairs_all_duplicates
echo pairs_all_duplicates_p=$pairs_all_duplicates_p
echo pairs_opt_duplicates=$pairs_opt_duplicates
echo pairs_opt_duplicates_p=$pairs_opt_duplicates_p
echo pairs_PCR_duplicates=$pairs_PCR_duplicates
echo pairs_PCR_duplicates_p=$pairs_PCR_duplicates_p

# Calculate Hi-C metrics
pairs_intra=`grep "number_of_intrachromosomal_pairs" $feather_qc | awk '{print $2}'`
pairs_intra_p=`echo "scale=4; 100 * $pairs_intra / $pairs_mapped_deduped" | bc | awk '{ printf("%.1f", $0) }'`
pairs_inter=`grep "number_of_interchromosomal_pairs" $feather_qc | awk '{print $2}'`
pairs_inter_p=`echo "scale=4; 100 * $pairs_inter / $pairs_mapped_deduped" | bc | awk '{ printf("%.1f", $0) }'`

echo pairs_intra=$pairs_intra
echo pairs_intra_p=$pairs_intra_p
echo pairs_inter=$pairs_inter
echo pairs_inter_p=$pairs_inter_p

pairs_intra_ge_15kb=`samtools view $bam_file_sorted | awk '{ if($7=="=" && $9 >= 15000) intra_ge_15kb++ } END {print intra_ge_15kb}'`
pairs_intra_ge_15kb_p=`echo "scale=4; 100 * $pairs_intra_ge_15kb / $pairs_mapped_deduped" | bc | awk '{ printf("%.1f", $0) }'`

echo pairs_intra_ge_15kb=$pairs_intra_ge_15kb
echo pairs_intra_ge_15kb_p=$pairs_intra_ge_15kb_p

# ChIP statistics
total_peaks=`wc -l $macs2_filepath | awk '{print $1}'`
bins_with_peaks=`awk '{ bin1 = int($2 / 5000); bin2 = int($3 / 5000); for(i = bin1; i <= bin2; i++) print $1"_"i }' $macs2_filepath | sort | uniq | wc -l`

short_VIPs_file=$feather_output"/"$dataset_name".shrt.vip.bed"
short_VIPs_file_sorted=$feather_output"/"$dataset_name".shrt.vip.sort.bed"

### Added support for multiple samples
cat $feather_output"/"*"."*".shrt.vip.bed" > $feather_output"/tmp.bed"
mv $feather_output"/tmp.bed" $short_VIPs_file

short_VIPs=`wc -l $short_VIPs_file | awk '{print $1}'`
cut -f1-3 $short_VIPs_file | sort -k1,1 -k2,2n -k3,3n > $short_VIPs_file_sorted
short_VIPs_with_peaks=`bedtools intersect -a $short_VIPs_file_sorted -b $macs2_filepath | wc -l`
short_VIPs_with_peaks_p=`echo "scale=4; 100 * $short_VIPs_with_peaks / $short_VIPs" | bc | awk '{ printf("%.1f", $0) }'`

echo total_peaks=$total_peaks
echo bins_with_peaks=$bins_with_peaks
echo short_VIPs=$short_VIPs
echo short_VIPs_with_peaks=$short_VIPs_with_peaks
echo short_VIPs_with_peaks_p=$short_VIPs_with_peaks_p

# HiChIP statistics
AND=`grep AND $maps_qc | awk '{sum=sum+$2}; END {print sum}'`
XOR=`grep XOR $maps_qc | awk '{sum=sum+$2}; END {print sum}'`
NOT=`grep NOT $maps_qc | awk '{sum=sum+$2}; END {print sum}'`

# New FRIPs calculation
long_intra_file=$feather_output"/"$dataset_name".long.intra.bedpe"
long_intra_filtered_file=$feather_output"/"$dataset_name".long.intra.filtered.bedpe"

### Added support for multiple samples
cat $feather_output"/"*"."*".long.intra.bedpe" > $feather_output"/tmp.bedpe"
mv $feather_output"/tmp.bedpe" $long_intra_file

cat $long_intra_file | awk -v resolution=$resolution 'BEGIN {resolution=resolution * 1000 } {R1_midpoint=(($2+$3)/2); R1_bin=int(R1_midpoint / resolution); R1_bin_start=(R1_bin * resolution); R1_bin_end=((R1_bin * resolution) + resolution - 1); R1_bin_mid=((R1_bin_start + R1_bin_end)/2); R2_midpoint=(($5+$6)/2); R2_bin=int(R2_midpoint / resolution); R2_bin_start=(R2_bin * resolution); R2_bin_end=((R2_bin * resolution) + resolution - 1); R2_bin_mid=((R2_bin_start + R2_bin_end)/2); distance=(sqrt((R1_bin_mid - R2_bin_mid)^2)); {if(distance >= 10000 && distance <= 2000000) { print $1"\t"R1_bin_start"\t"R1_bin_end"\t"$4"\t"R2_bin_start"\t"R2_bin_end } } }' | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n > $long_intra_filtered_file

# Calculate statistics based on Arima's AND, XOR and NOT metrics as well as target raw PE depth
ALL_Arima=`cat $long_intra_filtered_file | wc -l`
AND_Arima=`cut -f4,15 $maps_output"/reg_raw.chr"*"."*"."$resolution"k.and" | grep -v "dist" | awk '{ if($2 > 1) a+=$1 } END {print a}'`
XOR_Arima=`cut -f4,15 $maps_output"/reg_raw.chr"*"."*"."$resolution"k.xor" | grep -v "dist" | awk '{ if($2 > 1) a+=$1 } END {print a}'`
FRIPs=$(( $AND_Arima + $XOR_Arima ))
NOT_Arima=$(( $ALL_Arima - $FRIPs ))
FRIPs_p=`echo "scale=4; 100 * $FRIPs / $ALL_Arima" | bc | awk '{ printf("%.1f", $0) }'`

FRIPs_per_peak_bins=`echo "scale=4; $FRIPs / $bins_with_peaks" | bc | awk '{ printf("%.1f", $0) }'`
FRIPs_per_peak_bins_per_million_deduped=`echo "scale=5; ($FRIPs / $bins_with_peaks) * 1000000 / $pairs_mapped_deduped" | bc | awk '{ printf("%.2f", $0) }'`
target_FRIPs_per_peak_bin=382
mapping_rate="0.85" #v1.9 update
sequence_uniqueness_rate="0.44" #v1.9 update
rounding_level="10000000"
target_raw_pairs=`echo "scale=6; ( $target_FRIPs_per_peak_bin / (0.0000007568 * $FRIPs_per_peak_bins_per_million_deduped) / ($mapping_rate * $sequence_uniqueness_rate) )" | bc | awk -v rounding_level=$rounding_level '{rounded_num = ((int($0 / rounding_level) + 1) * rounding_level)} { printf("%d", rounded_num) }'` #v1.9 update

echo AND=$AND
echo XOR=$XOR
echo NOT=$NOT
echo ALL_Arima=$ALL_Arima
echo AND_Arima=$AND_Arima
echo XOR_Arima=$XOR_Arima
echo NOT_Arima=$NOT_Arima
echo FRIPs=$FRIPs
echo FRIPs_p=$FRIPs_p
echo FRIPs_per_peak_bins=$FRIPs_per_peak_bins
echo FRIPs_per_peak_bins_per_million_deduped=$FRIPs_per_peak_bins_per_million_deduped
echo target_raw_pairs=$target_raw_pairs

# Write QC tables

QC_result_deep=$outdir"/"$dataset_name"_Arima_QC_deep.txt"
QC_result_shallow=$outdir"/"$dataset_name"_Arima_QC_shallow.txt"

header_deep=("Sample Name" "Raw PE Reads" "Mapped Reads" "% Mapped Reads" "PCR Dups" "% PCR Dups" "Optical Dups" "% Optical Dups" "Mapped and De-duped" "% Mapped and De-duped" "Loops" "LR FRIPS per Peak Bin per Million Deduped PE Reads" "LR FRIPS per Peak Bin" "INTRA pairs" "% INTRA pairs" "INTRA > 15kb pairs" "% INTRA > 15kb pairs" "INTER pairs" "% INTER pairs" "ChIP Peaks" "Peak Bins" "Short VIPs" "Short VIPs in Peaks" "% Short VIPs in Peaks" "Enrichment")
IFS=$'\t'; echo "${header_deep[*]}" > $QC_result_deep

result_deep=($dataset_name $pairs_raw $pairs_uniq_mapped_MAPQ_ge_30 $pairs_uniq_mapped_MAPQ_ge_30_p $pairs_PCR_duplicates $pairs_PCR_duplicates_p $pairs_opt_duplicates $pairs_opt_duplicates_p $pairs_mapped_deduped $pairs_mapped_deduped_p $total_loops $FRIPs_per_peak_bins_per_million_deduped $FRIPs_per_peak_bins $pairs_intra $pairs_intra_p $pairs_intra_ge_15kb $pairs_intra_ge_15kb_p $pairs_inter $pairs_inter_p $total_peaks $bins_with_peaks $short_VIPs $short_VIPs_with_peaks $short_VIPs_with_peaks_p $enrichment)
IFS=$'\t'; echo "${result_deep[*]}" >> $QC_result_deep

header_shallow=("Sample Name" "Raw PE Reads" "Mapped Reads" "% Mapped Reads" "PCR Dups" "% PCR Dups" "Optical Dups" "% Optical Dups" "Mapped and De-duped" "% Mapped and De-duped" "Target Raw PE Reads" "LR FRIPS per Peak Bin per Million Deduped PE Reads" "LR FRIPS per Peak Bin" "INTRA pairs" "% INTRA pairs" "INTRA > 15kb pairs" "% INTRA > 15kb pairs" "INTER pairs" "% INTER pairs" "ChIP Peaks" "Peak Bins" "Short VIPs" "Short VIPs in Peaks" "% Short VIPs in Peaks" "Enrichment")
IFS=$'\t'; echo "${header_shallow[*]}" > $QC_result_shallow

result_shallow=($dataset_name $pairs_raw $pairs_uniq_mapped_MAPQ_ge_30 $pairs_uniq_mapped_MAPQ_ge_30_p $pairs_PCR_duplicates $pairs_PCR_duplicates_p $pairs_opt_duplicates $pairs_opt_duplicates_p $pairs_mapped_deduped $pairs_mapped_deduped_p $target_raw_pairs $FRIPs_per_peak_bins_per_million_deduped $FRIPs_per_peak_bins $pairs_intra $pairs_intra_p $pairs_intra_ge_15kb $pairs_intra_ge_15kb_p $pairs_inter $pairs_inter_p $total_peaks $bins_with_peaks $short_VIPs $short_VIPs_with_peaks $short_VIPs_with_peaks_p $enrichment)
IFS=$'\t'; echo "${result_shallow[*]}" >> $QC_result_shallow

echo -e "\nMAPS and Arima QC pipeline finished successfully!"
echo "Please download the QC result from: $QC_result_deep $QC_result_shallow and then copy the contents to the corresponding tables in the QC worksheet."
