![alt text](https://arimagenomics.com/public/images/header-logo.png "Celebrating Science and Scientist")


# MAPS Pipeline for Arima-HiChIP with Arima-HiC<sup>+</sup> Kit

This pipeline is for analyzing HiChIP data with MAPS and for generating Arima Genomics QC metrics and output data files.  The pipeline runs the normal MAPS algorithm from https://github.com/ijuric/MAPS with some minor changes to the default parameters. These parameters have been optimized by internal benchmarking and have been found to optimize sensitivity and specificity.  This pipeline will also generate shallow and deep sequencing QC metrics which can be copied into the Arima-HiChIP_QC_Worsheet.xls for analysis.  Additionally, the pipeline automatically generates metaplots for data QC and arc plots of chromatin loops for data visualization.

## Getting Started

This pipeline is suited for HiChIP from any protocol but optimal results are obtained when using the Arima-HiC<sup>+</sup> kit with the Arima-HiChIP protocol.

To order Arima-HiC<sup>+</sup> kits, please visit our website:

https://arimagenomics.com/

### Installing MAPS

Refer to the documentation at: https://github.com/ijuric/MAPS for help installing MAPS and its dependencies.

```
git clone https://github.com/ijuric/MAPS.git
```

The Arima genomic_features files and run_pipeline_Arima_v1.9.sh will be in the "MAPS/Arima_Genomics/" directory after downloading.

### Dependencies

#### python 3.4 (or later)

- Anaconda3

```
curl -O https://repo.anaconda.com/archive/Anaconda3-2019.07-Linux-x86_64.sh
chmod 777 Anaconda3-2019.07-Linux-x86_64.sh
./Anaconda3-2019.07-Linux-x86_64.sh
source ~/anaconda3/bin/activate
conda init
# read the user agreement and answer yes. install Anaconda3 in the default location which is your home folder.
```

- deeptools (v3.3.0), https://deeptools.readthedocs.io/en/develop/content/installation.html

```
conda install -c bioconda deeptools
```

- pandas (v0.20.3)

```
conda install pandas=0.20.3
```

- numpy (v1.13.1)

```
conda install numpy=1.13.1
```

- pysam (v0.15.2),

```
conda install pysam
```

- pybedtools (v0.8.0),

```
conda install --channel conda-forge --channel bioconda pybedtools
```

#### R 3.4.3

- argparse (v2.0.1)

```
# in R
install.packages("argparse")
```


#### Other

- bedtools (v2.27.1), https://bedtools.readthedocs.io/en/latest/content/installation.html

```
wget https://github.com/arq5x/bedtools2/archive/v2.27.1.tar.gz
tar -zxvf v2.27.1.tar.gz
cd bedtools2
make
```

- HTSLIB (v1.10.2), http://www.htslib.org/download/

```
wget https://github.com/samtools/htslib/releases/download/1.10.2/htslib-1.10.2.tar.bz2
tar -vxjf htslib-1.10.2.tar.bz2
cd ~/htslib-1.10.2/
./configure --prefix=/where/to/install
make
make install
export PATH=/where/to/install/bin:$PATH
```

- samtools (v1.10), http://www.htslib.org/download/

```
wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2
tar -vxjf samtools-1.10.tar.bz2
cd ~/samtools-1.10/
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

## Arima Specific Outputs

### Arima Shallow Sequencing QC

#### [output directory]/Arima_QC_shallow.txt
Contents: This file includes QC metrics for assessing the shallow sequencing data for each HiCHIP library.
- Break down of the number of read pairs
- The target sequencing depth for deep sequencing
- The breakdown of MAPS reads used for loop calling (AND and XOR reads)
- Summary statistics of the ChIP peak file used for loop calling
- The percentage of Short VIP's that overlap the ChIP peaks

### Arima Deep Sequencing QC

#### [output directory]/Arima_QC_deep.txt
Contents: This file includes QC metrics for assessing the deep sequencing data for each HiCHIP library.
- Break down of the number of read pairs
- The number of loops called
- The breakdown of MAPS reads used for loop calling (AND and XOR reads)
- Summary statistics of the ChIP peak file used for loop calling
- The percentage of Short VIP's that overlap the ChIP peaks

### Arima Arc Plots

#### [output directory]/arcplot_and_metaplot/[sample name].arcplot.gz
- bgzipped arc plot file

#### [output directory]/arcplot_and_metaplot/[sample name].arcplot.gz.tbi
- tabix index of the bgzipped arc plot file

These files can be viewed in the WashU EpiGenome Browser (http://epigenomegateway.wustl.edu/browser/).  See the Arima-HiChIP Analysis User Guide for more details.

### Metaplots, Meatmaps, and 1-D ChIP Signal

#### [output directory]/arcplot_and_metaplot/[sample name].bigwig
- bigwig file of all mapped reads.  This file can be used to view the sequencing coverage in a genome browser.

#### [output directory]/arcplot_and_metaplot/[sample name].heatmap.pdf
- PDF of Metaplot of all mapped reads overlapping the ChIP peaks used for loop calling.  This file can be used to assess the signal to noise of the HiChIP enrichment


## Arima Test Data

Test data to validate proper installation of the Arima MAPS pipeline can be found in [PATH_TO_MAPS]/MAPS/Arima_Genomics/test_data/ and includes:
- Test dataset (5M paired-reads) Arima-MAPS-test_R1.fastq.gz and Arima-MAPS-test_R1.fastq.gz
- ChIP peak file to be used to analysis the above data using the Arima-MAPS pipeline: ENCFF247YHM.UW.bed
- QC tables output by Arima-MAPS: Arima-MAPS-test_Arima_QC_shallow.txt and Arima-MAPS-test_Arima_QC_deep.txt as well as the same data in the Arima-HiChIP QC Worksheet that accompanies the Arima-HiChIP User Guide for Mammalian cell lines: Arima-HiChIP_QC_Worksheet_MAPS1.9_Test_data.xlsx
- Heatmap and Metaplot output by Arima-MAPS: Arima-MAPS-test.heatmap.pdf
- 1-D ChIP enrichment signal output by Arima-MAPS: Arima-MAPS-test.coverage.bigwig
- Loop calls from Arima-MAPS: Arima-MAPS-test.5k.2.sig3Dinteractions.bedpe
- Arc Plots generated from Arima-MAPS loop calls which can be viewed in the WashU Epigenome Browser (http://epigenomegateway.wustl.edu/browser/): Arima-MAPS-test.5k.2.arcplot.gz and Arima-MAPS-test.5k.2.arcplot.gz.tbi

The Arima-HiChIP Bioinformatics User Guide walks through an example of how to run the Arima-MAPS pipeline using this test data and provides additional information on the output files.

## Arima Pipeline Version

1.9

## Support

Please contact Ming Hu at hum@ccf.org for support using MAPS

For Arima customer support, please contact techsupport@arimagenomics.com

## Acknowledgments

Thank you to Ming Hu, Armen Abnousi, and Ivan Juric for collaborating with us on this pipeline and for hosting it on their GitHub page.
