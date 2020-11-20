![alt text](https://arimagenomics.com/public/images/header-logo.png "Celebrating Science and Scientist")


# MAPS Pipeline for Arima-HiChIP with Arima-HiC<sup>+</sup> Kit
This pipeline is for analyzing HiChIP data with MAPS and for generating Arima Genomics QC metrics and output data files. The pipeline runs the normal MAPS algorithm from https://github.com/ijuric/MAPS with some minor changes to the default parameters. These parameters have been optimized by internal benchmarking and have been found to optimize sensitivity and specificity. This pipeline will also generate shallow and deep sequencing QC metrics which can be copied into the Arima-HiChIP_QC_Worsheet.xls for analysis. Additionally, the pipeline automatically generates metaplots for data QC and arc plots of chromatin loops for data visualization.

## Getting Started
This pipeline is suited for HiChIP from any protocol but optimal results are obtained when using the Arima-HiC<sup>+</sup> kit with the Arima-HiChIP protocol.

To order Arima-HiC<sup>+</sup> kits, please visit our website:

https://arimagenomics.com/

### Installing MAPS
Refer to the documentation at: https://github.com/ijuric/MAPS for help installing MAPS and its dependencies.

```
git clone https://github.com/ijuric/MAPS.git
```

The Arima genomic_features files and run_pipeline_Arima_v2.0.sh will be in the "MAPS/Arima_Genomics/" directory after downloading.

### Dependencies

#### GCC Compiler (including gcc, g++ and make)

```
sudo apt update
sudo apt install build-essential
```

#### python 3.4 (or later)

- Anaconda3

```
curl -O https://repo.anaconda.com/archive/Anaconda3-2019.07-Linux-x86_64.sh
sh Anaconda3-2019.07-Linux-x86_64.sh
# read the user agreement and answer yes. install Anaconda3 in the default location which is your home folder.
```

- deeptools (v3.3.0 or later), https://deeptools.readthedocs.io/en/develop/content/installation.html

```
conda install -c bioconda deeptools
```

- pandas (v0.20.3 or later)

```
conda install pandas
```

- numpy (v1.13.1 or later)

```
conda install numpy
```

- pysam (v0.15.2 or later)

```
conda install pysam
```

- pybedtools (v0.8.0 or later)

```
conda install -c conda-forge -c bioconda pybedtools
```

- MACS2 (v2.2.7) https://github.com/macs3-project/MACS

```
conda install -c bioconda macs2
```

#### R 3.4.s or later

```
conda install R
```

- argparse (v2.0.1), https://cran.r-project.org/web/packages/argparse/index.html

```
# in R
install.packages("argparse")
```

- VGAM (v1.1-2), https://cran.r-project.org/web/packages/VGAM/index.html

```
# in R
install.packages("VGAM")
```

- data.table (v1.12.8), https://cran.r-project.org/web/packages/data.table/index.html

```
# in R
install.packages("data.table")
```

#### Other

- bedtools (v2.27.1), https://bedtools.readthedocs.io/en/latest/content/installation.html

```
conda install bedtools
```

- HTSLIB (v1.10.2), http://www.htslib.org/download/

```
wget https://github.com/samtools/htslib/releases/download/1.10.2/htslib-1.10.2.tar.bz2
tar -vxjf htslib-1.10.2.tar.bz2
cd htslib-1.10.2/
./configure --prefix=/where/to/install
make
make install
export PATH=/where/to/install/bin:$PATH
```

- samtools (v1.10), http://www.htslib.org/download/

```
wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2
tar -vxjf samtools-1.10.tar.bz2
cd samtools-1.10/
./configure --prefix=/where/to/install
make
make install
export PATH=/where/to/install/bin:$PATH
```

- bcftools (v1.10.2), http://www.htslib.org/download/

```
wget https://github.com/samtools/bcftools/releases/download/1.10.2/bcftools-1.10.2.tar.bz2
tar -vxjf bcftools-1.10.2.tar.bz2
cd bcftools-1.10.2
./configure --prefix=/where/to/install
make
make install
export PATH=/where/to/install/bin:$PATH
```

- bwa (v0.7.17) https://github.com/lh3/bwa

```
conda install -c bioconda bwa
```

#### Add python, R and MACS2 (optional) installation paths to $PATH variable

```
PATH=$PATH:/path_to_Python/:/path_to_R/:[:/path_to_MACS2/]
```

## Usage (Command line options)
Arima-MAPS_v2.0.sh [-C call_peaks] [-p peak_type] [-F feather] [-M maps]
       [-I fastq_dir_with_file_prefix] [-O outdir] [-m macs2_filepath]
       [-o organism] [-b bwa_index] [-t threads] [-f patterned_flowcell]
       [-P plot] [-H generate_hic] [-s bin_size] [-r binning_range] [-d fdr]
       [-Q mapq] [-l length_cutoff] [-h]

       * [-C call_peaks]: 0 (default) to use ChIP peak file provided, 1 to call peaks
           using MACS2
       * [-p peak_type]: broadness of chromatin factor. Required if call_peaks=1. Must
           be either "broad" (H3K4me3) or "narrow" (CTCF). Both choices result in broad
           peaks being called by MACS2.
       * [-F feather]: 1 (default) to run feather, 0 to skip
       * [-M maps]: 1 (default) to run MAPS on data processed with feather, 0 to skip
       * [-I fastq_dir_with_file_prefix]: absolute path and file prefix of the fastq
           files, up to "_R1" or ".R1"
       * [-O outdir]: directory which the feather and MAPS outputs will be placed
       * [-m macs2_filepath]: path to the ChIP peak file that will be used for loop
           calling with MAPS. Required if call_peaks=0. Ignore or leave it empty ("")
           if call_peaks=1
       * [-o organism]: organism of the genomic feature file to be used, options:
           "mm9", "mm10", "hg19", and "hg38"
       * [-b bwa_index]: absolute path to the reference genome sequence (.fa) which will
           be used to derive the BWA index as well. Ex: "/home/reference_sequence/hg19.fa"
       * [-t threads]: number of threads for MAPS to run on. Use 4-8 for shallow sequencing
           and 12-20 for deep sequencing.
       * [-f patterned_flowcell]: Use 1 for deep sequencing and 0 for shallow sequencing
           datasets. "1" if the data was sequenced on a patterned flowcell, "0" if
           non-patterned flowcell. This is used for calculating optical duplicates and
           PCR duplicate rates.
       * [-P plot]: 1 (default) to generate plots, 0 to skip
       * [-H generate_hic]: 1 (default) to generate a .hic file for visualization with
           Juicer, 0 to skip
       * [-s bin_size]: resolution of the loops called. Default: 5000
       * [-r binning_range]: Maximum distance for loop calling. Default: 2000000
       * [-d fdr]: Do not change! False Discovery Rate threshold: 1 = 0.1,
           2 (default) = 0.01, 3 = 0.001, ect...
       * [-Q mapq]: Phred scaled mapping quality threshold. Default: 30
       * [-l length_cutoff]: Minimum genomic distance. Default: 1000
       * [-h]: print this help and exit

### Example

- With public / known ChIP peaks

```
[PATH TO MAPS INSTALLATION]/bin/Arima-MAPS_v2.0.sh -C 0 -I [FASTQ DIR]/Arima-MAPS-test -O [PATH TO OUTPUT DIR] -m [PEAK FILE DIR]/ENCFF247YHM.UW.bed -o hg19 -b [GENOME DIR]/hg19.fa -t 8 -f 1
```

- With MACS2 peak calls

```
[PATH TO MAPS INSTALLATION]/bin/Arima-MAPS_v2.0.sh -C 1 -p broad -I [FASTQ DIR]/Arima-MAPS-test -O [PATH TO OUTPUT DIR] -o hg19 -b [GENOME DIR]/hg19.fa -t 8 -f 1
```

## Arima Specific Outputs

### Arima Shallow Sequencing QC

#### [output directory]/[sample name]_Arima_QC_shallow.txt
Contents: This file includes QC metrics for assessing the shallow sequencing data for each HiChIP library.
- Break down of the number of read pairs
- The target sequencing depth for deep sequencing
- The breakdown of MAPS reads used for loop calling (AND and XOR reads)
- Summary statistics of the ChIP peak file used for loop calling
- The percentage of Short VIP's that overlap the ChIP peaks
- The peak enrichment score

### Arima Deep Sequencing QC

#### [output directory]/[sample name]_Arima_QC_deep.txt
Contents: This file includes QC metrics for assessing the deep sequencing data for each HiChIP library.
- Break down of the number of read pairs
- The number of loops called
- The breakdown of MAPS reads used for loop calling (AND and XOR reads)
- Summary statistics of the ChIP peak file used for loop calling
- The percentage of Short VIP's that overlap the ChIP peaks
- The peak enrichment score

### Arima Arc Plots

#### [output directory]/arcplot_and_metaplot/[sample name].*.arcplot.gz
- bgzipped arc plot file

#### [output directory]/arcplot_and_metaplot/[sample name].*.arcplot.gz.tbi
- tabix index of the bgzipped arc plot file

These files can be viewed in the WashU EpiGenome Browser (http://epigenomegateway.wustl.edu/browser/). See the Arima-HiChIP Analysis User Guide for more details.

### Metaplots, heatmaps, and 1-D ChIP Signal

#### [output directory]/arcplot_and_metaplot/[sample name].coverage.bigwig
- bigwig file of all mapped reads. This file can be used to view the sequencing coverage in a genome browser.

#### [output directory]/arcplot_and_metaplot/[sample name].heatmap.pdf
- PDF of Metaplot of all mapped reads overlapping the ChIP peaks used for loop calling. This file can be used to assess the signal to noise of the HiChIP enrichment.

### Peaks generated using MACS2 (only available when call_peaks=1)

#### [output directory]/MACS2_peaks/[sample name]_peaks.broadPeak
- ChIP peak file generated using MACS2

#### [output directory]/MACS2_peaks/[sample name]_peaks.xls
- ChIP peak file generated using MACS2 in Microsoft Excel format

#### [output directory]/MACS2_peaks/[sample name]_model.r
- An R script for generating the peak model

## Arima Test Data

Test data to validate proper installation of the Arima MAPS pipeline can downloaded from our ftp site by running the following command:

```
wget ftp://ftp-arimagenomics.sdsc.edu/pub/MAPS/test_data/*
```

The Download includes:
- Test dataset (5M paired-reads) Arima-MAPS-test_R1.fastq.gz and Arima-MAPS-test_R2.fastq.gz
- ChIP peak file to be used to analyze the above data using the Arima-MAPS pipeline: ENCFF247YHM.UW.bed
- QC tables output by Arima-MAPS: Arima-MAPS-test_Arima_QC_shallow.txt and Arima-MAPS-test_Arima_QC_deep.txt as well as the same data in the Arima-HiChIP QC Worksheet that accompanies the Arima-HiChIP User Guide for Mammalian cell lines: Arima-HiChIP_QC_Worksheet_MAPS2.0_Test_data.xlsx
- Heatmap and Metaplot output by Arima-MAPS: Arima-MAPS-test.heatmap.pdf
- 1-D ChIP enrichment signal output by Arima-MAPS: Arima-MAPS-test.coverage.bigwig
- Loop calls from Arima-MAPS: Arima-MAPS-test.5k.2.sig3Dinteractions.bedpe
- Arc Plots generated from Arima-MAPS loop calls which can be viewed in the WashU Epigenome Browser (http://epigenomegateway.wustl.edu/browser/): Arima-MAPS-test.5k.2.arcplot.gz and Arima-MAPS-test.5k.2.arcplot.gz.tbi

The Arima-HiChIP Bioinformatics User Guide walks through an example of how to run the Arima-MAPS pipeline using this test data and provides additional information on the output files.

## Arima Pipeline Version

2.0

## Whats New in this Version

- Added command line options.
- Added enrichment score to the QC table.

## Support

Please contact Ming Hu at hum@ccf.org for support using MAPS

For Arima customer support, please contact techsupport@arimagenomics.com

## Acknowledgments

Thank you to Ming Hu, Armen Abnousi, and Ivan Juric for collaborating with us on this pipeline and for hosting it on their GitHub page.
