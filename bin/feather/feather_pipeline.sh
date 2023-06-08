#run using python3 (tried on python 3.6)
#pysam matching with samtools 1.7+
#accepts .sam/.bam/.fastq or .fastq.gz for fastq1 and fastq2 variables below

fastq1=
fastq2=
bwa_index=
outdir=
prefix=
mapq=
length_cutoff=
threads=8
memory_per_thread="20G"
per_chr=1 # set this to one if you don't want per chromosome output bed and bedpe files

./feather_pipe preprocess -o $outdir -p $prefix -f1 $fastq1 -f2 $fastq2 -b $bwa_index -q $mapq -l $length_cutoff -t $threads -m $memory_per_thread -c $per_chr
