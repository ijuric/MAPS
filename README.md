# MAPS pipeline user manual
## What is MAPS pipeline?
MAPS (Model-based Analysis of PLAC-Seq data) pipeline is a a set of multiple scripts used to analyze PLAC-Seq and HiChIP data. MAPS pipeline contains multiple scripts. The run_pipeline.sh is a prototypes of a shell script containing the locations of files needed to run MAPS is entered as well as calls for different MAPS pipeline scripts. Internally, there are two parts of MAPS pipeline. First part, called feather, does mapping and preprocessing of pair-end reads and creates long and short .bed/.bedpe files. Second part, called MAPS, does read binning and peak calling.
For additional description of MAPS method, check our paper: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006982
For any questions regarding MAPS please send mail to Ivan Juric ( ivan.juric.gen@gmail.com ) or Ming Hu ( hum@ccf.org )

**MAPS runs with Arima HiChIP kit using a slightly modified run_pipeline file. Please refer to https://github.com/ijuric/MAPS/tree/master/Arima_Genomics for more information**

## MAPS pipeline setup
### Downloading prerequisites and MAPS pipeline
MAPS requires following programs and packages. Install them prior to using MAPS. MAPS runs on Linux.

**python 3.4 (or later)**

Python libraries:
  * pandas (v0.20.3)*
  * numpy (v1.13.1)
  * itertools (v3.2)
  * pysam(v0.15.2)
  * pybedtools(v0.8.0)
  
**R 3.4.3**

R Packages:
  * MASS (v7.3-50)
  * VGAM package (v1.0.5)
  * data.table package (v1.11.2)
  * bedtools (v2.27.1)
  * Samtools(v1.7 or later)
  * Bwa (v0.7.12)

**Juicer tools (needed only if you want .hic file; in feather dir on github)**
https://github.com/theaidenlab/juicer/wiki/Download

*Numbers in parentheses are versions we used, and while it is very likely that MAPS will run with newer versions, we do not guarantee it.*

Add bwa, bedtools, samtools to the path. Check here how to add to set up path on linux:
http://hmgaudecker.github.io/econ-python-environment/paths.html

Index your reference genome using bwa
To index genome using bwa, type:
'''
bwa index [ref_genome_fastq_file]
'''
Download reference genome. In our paper, we use ones downloaded from ucsc genome browser. For example:
http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz

Alternatively, you can download genomes mentioned in Heng Li's blog here: https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use

MAPS requires files with genome features information to run properly. We provide those files for 5kb and 10kb resolution. If you want to run MAPS at different resolution, download appropriate genome features files from here:
http://enhancer.sdsc.edu/yunjiang/resources/genomic_features/

Download or clone from here: https://github.com/ijuric/MAPS

## Running MAPS
Run one run_pipeline script per biological replica. For each biological replica, run_pipeline script runs feather, MAPS and MAPS-tools parts. Run MAPS pipeline for each replica before running it for merged dataset. That way, we will not need to map reads for merged dataset, but will use .bed/.bedpe files from each replica. It is crucial to add appropriate information about all biological replicas in the run_pipeline script. We show below how to do that.

Copy MAPS/bin/run_pipeline.sh to MAPS/bin/run_pipeline_[PROJECT_NAME].sh, where [PROJECT_NAME] is the name of your set of biological replicas.

In run_pipeline_[PROJECT_NAME].sh set parameters:

* python_path
  * location of python (MAPS works with python3)
* Rscript_path
  * location of Rscript
* feather
  * 1 if you want to run feather part. 0 if not. You want this set to 0 if you have already ran feather for this dataset and need to run MAPS only. Otherwise set it to 1. If you are running MAPS on merged datasets, feather must be set to 1
* maps
  * 1 if you want to run MAPS past. 0 if not. By setting this to 0, you ran mapping and prefiltering (outputs are .bedpe and .bed files), but not binning or peak calling.
* number_of_datasets
  * Number of biological replicas. Can be 1,2, … N. If set to 1, script assumes that user is running MAPS for a particular biological replica. If number_of_datasets is larger than 1, script assumes that user is running MAPS on merged data, and it will look for .bedpe and .bed files of each replica defined in dataset1, dataset2,..., datasetN lines (See below)
* dataset_name
  * Name of dataset.
* fastq_format
  * Extension of .fastq files (usually .fastq, but can be .fa) Your fastq file names must be named [datast_name]_R1.fastq (or .fa) and [datast_name]_R2.fastq, where [dataset_name] is the name of your data set
* fastq_dir
  * Location of fastq files
* outdir
  * Location of MAPS output files
* macs2_filepath
  * Location of MACS2 peaks
* organism
  * either mm10, hg19 or hg38, depending on which chromosome you use. Important for genomic features file. Also selecting this defines the number of chromosomes (19 for mouse, 22 for human).
* bwa_index
  * index of bwa. Set to path to bwa indexed genome
* bin_size
  * resolution. Usually 5000 or 10000.  
* binning_range
  * binning range. How far 3D interactions can be called, also affects the estimate of the expected count. Default=1000000 Do not set to high value if data is sparse. Check MAPS paper for more details
* filter_file
  * Location of file with blacklisted bins. This is used if you want to exclude genomic regions from MAPS analysis. Reads mapping to those regions will be ignored. Set to “None” if not blacklisting anything. A blacklist file is a tab delimited table containing two columns, named chr and bin, representing chromosome and the start location of bins that user wants to blacklist. For example:
```
chr	bin
chr1	10000
chr1	15000
chr13	7000000
chr15	1300000
```
* generate_hic
  * 1 if you want .hic file to be generated, 0 if not.
* Dataset1, dataset2,...
  * Used for running MAPS on merged dataset. When number_of_datasets > 1, run_pipeline looks at the paths defined by Dataset1, dataset2… to find  where each biological replica .bedpe and .bed files are.
* model
  * regression model. Can be pospoisson or negbinom. MAPS supports positive poisson (pospoisson) or negative binomial (negbinom) regression. Default value is pospoisson (positive poisson regression).
* sex_chroms_to_process
  * either X,Y,XY or NA. This specifies which (if any) sex chromosomes the user wants to run MAPS on: X = X chr only, Y = Y chr only, XY = both X and Y chroms, NA = none (just autosomal). Default value is NA.
* bwa_index
  * path to bwa indexed reference genome
* fdr
  * this is used for labeling purposes only. Do not change. See the MAPS paper for more details on how fdr is picked.

### Running MAPS
Run appropriately set up run_pipeline script
## MAPS Output
MAPS generates multiple files. Main output is .peaks.bedpe file which contains the list of significant 3D interactions in .bedpe format.

* chr1
  * bin1 chromosome
* start1
  * bin1 start position
* end1
  * bin1 end position
* chr2
  * bin2 chromosome
* start2
  * bin2 start position
* end2
  * bin2 end position
* count
  * the number of pair-end reads observed for [bin1, bin2] pair (raw count)
* expected
  * the expected count for [bin1, bin2] pair (regression result)
* fdr
  * false discovery rate associated with count/expected ratio
* ClusterLabel
  * the name (label) of the cluster to which this bin pair belongs to
* ClusterSize
  * the size of the cluster to which this bin pair belongs to
* ClusterType
  * the type of cluster (Singleton, Sharp Peak, Broad Peak) of cluster to which this bin pair belongs to.
* ClusterNegLog10P
  * -log10(P-value) of the cluster to which this bin pair belongs to
* ClusterSummit
  * 1 is this bin pair is cluster summit, 0 otherwise. A bin pair b[i,j] is cluster summit if no other bins in this cluster have lower fdr value than b[i,j]. Note that cluster can have multiple summits if multiple bin pairs have the same (lowest) fdr value.


## Example scripts ##

### test example

File bin/run_pipeline_test.sh is a small examlpe that runs fast. Use it to confirm that MAPS is running properly. Currently, it is set to run from my home directory (/home/jurici/). You will need to change the following paths to correspond to locations of files on your computer:

```
python_path=/home/abnousa/software/python3.6.5/bin/python #should have pysam, pybedtools installed. bedtools, samtools should be in the path
Rscript_path=/opt/R-3.4.3/lib64/R/bin/Rscript
fastq_dir="/home/jurici/MAPS/examples/test_set1"
outdir="/home/jurici/MAPS/examples/test_set1/output"
macs2_filepath="/home/jurici/MAPS/examples/test_set1/macs2_peaks_final.replicated.narrowPeak"
bwa_index="/home/jurici/MAPS/MAPS_data_files/"$organism"/BWA_index/mm10_chrAll.fa"
```
fastq and 1D peaks (macs2) files are in your MAPS folder (/home/[USER]/MAPS/examples)

Run run_pipeline_test.sh script

Results should looks omething like this:
https://github.com/ijuric/MAPS/blob/master/examples/expected_output_test_set1.bedpe

### GM12878 cells example:
In this example (bin/run_pipeline_example.sh), we show how to run MAPS on one HiChIP dataset on Smc1a cohesin subunit in human B lymphocyte GM12878 cells. After successfully finishing setting described in MAPS pipeline setup, you will need to download data files:
* 1D chip-seq peaks
  * https://www.encodeproject.org/files/ENCFF686FLD/@@download/ENCFF686FLD.bed.gz
* GSM2138324: HiChIP GM biological replicate 1 technical replicate 1 and technical replica 2 data:
  * https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR3467175
  * https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR3467176

To download HiChIP GM data from NCBI site, you will need NCBI SRA Toolkit, which can be downloaded here: https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software Once you download and set up SRA Toolkit, to download .fastq files type:
```
./fastq-dump --split-files SRR3467175
./fastq-dump --split-files SRR3467176
```
This will download SRR3467175_1.fastq, SRR3467175_2.fastq, SRR3467176_1.fastq and SRR3467176_2.fastq files to your computer. Next, merge fastq files from technical replicas: 
```
cat SRR3467175_1.fastq SRR3467176_1.fastq > GM12878_R1.fastq
cat SRR3467175_2.fastq SRR3467176_2.fastq > GM12878_R2.fastq
```
This is a proper way of merging fastq files from technical replicas becasue all duplicated reads will be removed during feather step of MAPS pipeline. We cannot merge biological replicas this way, since then we risk loosing some reads when removing duplicated. This is why, when we run MAPS on merged biologal replicas, we first need to run MAPS on each replica first, and then use mapped and filtered reads to call 3D interactions.

run run_pipeline_example.sh

After you run this example, output should look like this:
https://github.com/ijuric/MAPS/blob/master/examples/GM12878.5k.2.peaks.bedpe
