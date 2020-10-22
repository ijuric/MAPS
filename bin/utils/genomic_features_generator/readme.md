## Creating Genomic Features Files:

Genomic features files are required for running MAPS. We have provided pre-generated files for commonly used genomes: mm10, hg19, and hg38 (GRCh38 and for various resolutions and for MboI restriction enzyme. However, if you are working with other genomes/resolutions/enzymes you will need to generate new files. The scripts provided here will help you to do so.

###Software Prerequisites:
In order to use the utility provided here you will need the following:
1. [Snakemake](https://snakemake.readthedocs.io): Can be installed using pip
2. [bedtools](https://bedtools.readthedocs.io)
3. bigWigAverageOverBed: Can be downloaded from UCSC Genome Browsers "utilites download" page [here](http://hgdownload.soe.ucsc.edu/admin/exe/). Find the directory that matches your operating system (e.g. linux.x86_64/) and then locate the bigWigAverageOverBed file and download it. After download make sure it is runnable by running: `chmod +x bigWigAverageOverBed`.
If the bedtools and bigWigAverageOverBed are not in your system's path, you should either add them to the path or speficy the complete path in the *config.json* file. Similarly, if the snakemake is not in the path, you should either add it to the path or specify the complete path in the last step mentioned [below](#run_script) when you run the script.

### Data Requirements:
To generate a genomic features file, you should have three files ready:
- genome sequence file (.fa/.fastq file, each chromosome as a separate sequence in the file). These files usually can be downloaded from UCSC genome browser download page or from the encode project.
- A file containing chromosome sizes for the genome. These files also usually can be obtained from the same place you download the genome file.
- A *bigwig* file containing the mappability information for the genome. For some genomes these files can be downloaded from sources such as Encode. For the less common/newer genomes, you might need to generate this file yourself by some other software. We have used *gem-mappability* from the [GEM  library] (https://sourceforge.net/projects/gemlibrary/files/gem-library/) with read length set to 50.  The output file should be converted to bigwig format.

### How to run:
1. <a name="step_json"></a>Add the enzyme information in the *config.json*  file (if it is not already there):
  - In the *"re*" section add a new category, giving it a unique name. (the enzyme name?)
  - Specify sequence and cleavage position of the enzyme <sup></sup>[want to know more details on this?](#config_json)</sup></a>.
  
  Example: For HhaI (GCG'C) your entry will look like:
   ```
"HhaI": {
                        "site": "GCGC",
                        "pos": "3"
                }   
   ```
2. Modify lines 5-11 of the *Snakefile_multienzyme* file:
   ```
   genome="hg19" #specify name of the genome
   enzyme="MboI" #this is the name you used in step1 in the config.json file 
   bin_size="5Kb" #binsize can be an integer or an integer followed by a Kb or Mb suffix
   genome_fastq="/path/to/genome/sequence/file" #this is file #1 mentioned above
   chrom_sizes="/path/to/chromosome/sizes/file" #this is file #2 mentioned above
   mappability_input="/path/to/mappability/bigwig/file" #this is file #3 mentioned above
   ```
3.  <a name="run_script"></a>Run the scripts using the command below:
   ```
   snakemake -s Snakefile_multienzyme
   ``` 
   </li>

After running these three steps, for the given *genome*, *enzyme*, *bin_size* and the *outdir* that you set, the output genomic features file will be stored in: `<outdir>/<genome>/<enzyme>/F_GC_M_<enzyme>_<bin_size>_el.txt`

<a name="config_json">[More on requirements for adding restriction enzymes to the config.json file](#step_json)</a>: If you are using multiple enzymes, they should all be specified here separated by a comma(*,*). Similarly, cleavage locations for all of the enzymes should be specifed separated by a comma. Avoid *spaces* before and after the commas. You can also use python regular expressions to specify mutliple restriction enzymes. Check how we have added Arime enzymes in the *config.json* file. We have one entry for MboI enzyme under *arima* entry, separated by a comma from another entry. The second entry is a regular expression equivalent to *GANTC* (yes, *dot* in python regular expressions means any character). 
To generate features files for MNase-based experiments, you can set the *re* parameter in the *Snakefile_multienzyme* file with *Mnase*. If you are generating features files for multiple resolutions, start with the highest resolution (for example first run 5Kb, then 10Kb, and so on). 
### Credits:
This [set of] scripts is provided by [**Yunjiang Qiu**](https://scholar.google.com/citations?user=0IzF8KEAAAAJ&hl=en). (Hey Yunjiang! you rock :thumbsup:).  We have slightly modified it to allow generating files for experiments that involve multiple enzymes such as the one performed by Arima HiChIP kits.

