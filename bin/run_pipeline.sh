#!/bin/bash
python_path=/home/abnousa/software/python3.6.5/bin/python #should have pysam, pybedtools installed. bedtools, samtools should be in the path
Rscript_path=/opt/R-3.5.0/bin/Rscript
###################################################################
feather=1 #start from feather or run only MAPS
maps=1
number_of_datasets=1
dataset_name="test_set1"
fastq_format=".fastq"
fastq_dir="/home/jurici/MAPS/examples/test_set1"
outdir="/home/jurici/MAPS/examples/outputs/test_set1"
macs2_filepath="/home/jurici/MAPS/PLAC-Seq_datasets/test_dataset2/MACS2_peaks/final.replicated.narrowPeak"
organism="mm10"
bwa_index=""
bin_size=10000
binning_range=1000000
fdr=2 # this is used just for labeling. do not change
filter_file="None"
generate_hic=1
mapq=30
length_cutoff=1000
threads=4
model="pospoisson" #"negbinom"
sex_chroms_to_process="X"
optical_duplicate_distance=0
####################################################################
###SET THE VARIABLES AT THIS PORTION ONLY IF 
### number_of_datasets > 1 (merging exisitng datasets)
### specify as many datasets as required
####################################################################
dataset1="/home/jurici/MAPS/PLAC-Seq_datasets/test_dataset2/feather_output/test_set1_current"
dataset2="/home/jurici/MAPS/PLAC-Seq_datasets/test_dataset2/feather_output/test_set2_current"
dataset3=""
dataset4=""
#...
##################################################################
###SET THESE VARIABLES ONLY IF FEATHER = 0 AND YOU WANT TO RUN
###USING A SPECIFIC FEATHER OUTPUT RATHER THAN $datasetname_Current
###################################################################
feather_output_symlink=""
##################################################################

DATE=`date '+%Y%m%d_%H%M%S'`
#####Armen:
fastq1=$fastq_dir/$dataset_name"_R1"$fastq_format
fastq2=$fastq_dir/$dataset_name"_R2"$fastq_format
feather_output=$outdir"/feather_output/"$dataset_name"_"$DATE
if [ "$feather_output_symlink" == "" ]; then
	feather_output_symlink=$outdir"/feather_output/"$dataset_name"_current"
fi
resolution=$(bc <<< "$bin_size/1000")
per_chr='True' # set this to zero if you don't want per chromosome output bed and bedpe files
feather_logfile=$feather_output"/"$dataset_name".feather.log"
resolution=$(bc <<< "$bin_size/1000")
cwd="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"
hic_dir="tempfiles/hic_tempfiles"
if [ $organism == "mm10" ]; then
	if [ -z $bwa_index ]; then
        	bwa_index="/home/jurici/MAPS/MAPS_data_files/"$organism"/BWA_index/mm10_chrAll.fa"
	fi
	genomic_feat_filepath=$cwd"/../MAPS_data_files/"$organism"/genomic_features/F_GC_M_MboI_"$resolution"Kb_el.mm10.txt"
	chr_count=19
elif [ $organism == "mm9" ]; then
	if [ -z $bwa_index ]; then
		bwa_index="/home/jurici/MAPS/MAPS_data_files/"$organism"/BWA_index/mm9.fa"
	fi
	genomic_feat_filepath=$cwd"/../MAPS_data_files/"$organism"/genomic_features/F_GC_M_MboI_"$resolution"Kb_el.mm9.txt"
	chr_count=19
elif [ $organism == "hg19" ]; then
	if [ -z $bwa_index ]; then
		bwa_index="/home/jurici/MAPS/MAPS_data_files/"$organism"/BWA_index/hg19.fa"
	fi
	genomic_feat_filepath=$cwd"/../MAPS_data_files/"$organism"/genomic_features/F_GC_M_MboI_"$resolution"Kb_el.hg19.txt"
	chr_count=22
elif [ $organism == "hg38" ]; then
	if [ -z $bwa_index ]; then
		bwa_index="/home/jurici/MAPS/MAPS_data_files/"$organism"/BWA_index/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
	fi
	genomic_feat_filepath=$cwd"/../MAPS_data_files/"$organism"/genomic_features/F_GC_M_MboI_"$resolution"Kb_el.GRCh38.txt"
	chr_count=22
fi

####Ivan:"
if [[ $sex_chroms_to_process != "X" && $sex_chroms_to_process != "Y" && $sex_chroms_to_process != "XY" ]]; then
	sex_chroms_to_process="NA"
	sex_chroms=""
else
	sex_chroms=$sex_chroms_to_process
fi
long_bedpe_dir=$feather_output_symlink"/"
short_bed_dir=$feather_output_symlink"/"
maps_output=$outdir"/MAPS_output/"$dataset_name"_"$DATE"/"
maps_output_symlink=$outdir"/MAPS_output/"$dataset_name"_current"
#genomic_feat_filepath="/home/jurici/MAPS/MAPS_data_files/"$organism"/genomic_features/"$genomic_features_filename

if [ $feather -eq 1 ]; then
	mkdir -p $feather_output
	if [ $number_of_datasets -ge 2 ]; then
		ds_array=()
		hic_array=()
		for i in `seq $number_of_datasets`
		do
			ds_name="dataset$i"
			hic_name="$ds_name/$hic_dir"
			eval ds=\$$ds_name
			eval hic=\$$hic_name
			ds_array+=($ds)
			hic_array+=($hic)
		done
		#printf '%s\n' "${qc_array[@]}"
		$python_path $cwd/feather/combine_feathers.py -o $feather_output -p $dataset_name -a $macs2_filepath -d "${ds_array[@]}"
		if [ $generate_hic -eq 1 ]; then
			echo "$feather_output/$hic_dir"
			mkdir -p $feather_output"/"$hic_dir
			$cwd/feather/concat_hic.sh  $feather_output $dataset_name $hic_dir "${hic_array[@]}"
		fi
	else
		$python_path $cwd/feather/feather_pipe preprocess -o $feather_output -p $dataset_name -f1 $fastq1 -f2 $fastq2 -b $bwa_index -q $mapq -l $length_cutoff -t $threads -c $per_chr -j $generate_hic -a $macs2_filepath -d $optical_duplicate_distance
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
fi
	
if [ $maps -eq 1 ]; then
	mkdir -p $maps_output
	echo "$dataset_name $maps_output $macs2_filepath $genomic_feat_filepath $long_bedpe_dir $short_bed_dir $bin_size $chr_count $maps_output"
	$python_path $cwd/MAPS/make_maps_runfile.py $dataset_name $maps_output $macs2_filepath $genomic_feat_filepath $long_bedpe_dir $short_bed_dir $bin_size $chr_count $maps_output $sex_chroms_to_process --BINNING_RANGE $binning_range
	echo "first"
	$python_path $cwd/MAPS/MAPS.py $maps_output"maps_"$dataset_name".maps"
	echo "second"
	$Rscript_path $cwd/MAPS/MAPS_regression_and_peak_caller.r $maps_output $dataset_name"."$resolution"k" $bin_size $chr_count$sex_chroms $filter_file $model 
	$Rscript_path $cwd/MAPS/MAPS_peak_formatting.r $maps_output $dataset_name"."$resolution"k" $fdr $bin_size
	echo "third"
	cp "$(readlink -f $0)" $maps_output"/execution_script_copy"
	chmod 777 $maps_output
	ln -sfn $maps_output $maps_output_symlink
fi
