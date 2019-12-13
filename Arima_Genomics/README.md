# MAPS Pipeline for Arima-HiChIP with Arima-HiC<sup>+</sup> Kit

This pipeline is for analyzing HiChIP data with MAPS and for generating Arima Genomics QC metrics and output data files.  The pipeline runs the normal MAPS algorithm from https://github.com/ijuric/MAPS with some minor changes to the default parameters. These parameters have been optimized by internal benchmarking and have been found to optimize sensitivity and specificity.  This pipeline will also generate shallow and deep sequencing QC metrics which can be copied into the Arima-HiChIP_QC_Worsheet.xls for analysis.  Additionally, the pipeline automatically generates metaplots for data QC and arc plots of chromatin loops for data visualization.

## Getting Started

This pipeline is suited for HiChIP from any protocol but optimal results are obtained when using the Arima-HiC<sup>+</sup> kit with the Arima-HiChIP protocol.

To order Arima-HiC<sup>+</sup> kits, please visit our website:

https://arimagenomics.com/

### Installing MAPS

Refer to the documentation at: https://github.com/ijuric/MAPS for help installing MAPS and its dependencies

### Arima Specific Prerequisites

This pipeline uses all the same dependencies as MAPS as well as:

#### python 3.4 (or later)

- deeptools, https://deeptools.readthedocs.io/en/develop/content/installation.html

```
conda install -c bioconda deeptools
```

#### R 3.4.3

- argparse

```
install.packages("argparse")
```

#### Other

- picard 2.6.0, https://broadinstitute.github.io/picard/


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

### Metaplots

#### [output directory]/arcplot_and_metaplot/[sample name].bigwig
- bigwig file of all mapped reads.  This file can be used to view the sequencing coverage in a genome browser.

#### [output directory]/arcplot_and_metaplot/[sample name].metaplot.pdf
- PDF of Metaplot of all mapped reads overlapping the ChIP peaks used for loop calling.  This file can be used to assess the signal to noise of the HiChIP enrichment


## Arima Pipeline Version

1.8

## Support

Please contact Ming Hu at hum@ccf.org for support using MAPS

For Arima customer support, please contact techsupport@arimagenomics.com

## Acknowledgments

Thank you to Ming Hu, Armen Abnousi, and Ivan Juric for collaborating with us on this pipeline and for hosting it on their GitHub page.
