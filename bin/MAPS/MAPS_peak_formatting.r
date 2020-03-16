options("scipen"=999)
library(data.table)

args <- commandArgs(trailingOnly=TRUE)

INDIR = args[1]
SET = args[2]
FDR = args[3]
RESOLUTION = as.numeric(args[4])

#peaks_raw = read.table(args[1], header=T)
#peaks_raw = read.table('../results/mESC/MY_113.MY_115.5k.2.peaks',header=T, stringsAsFactors=F)

inf = paste(INDIR, SET,'.',FDR,'.peaks',sep='')
peaks_raw = read.table(inf,header=T, stringsAsFactors=F)

peaks = as.data.table(subset(peaks_raw, ClusterType != 'Singleton' | (ClusterType == 'Singleton' & fdr < 1e-4))) # remove singletons
peaks[, summit := 1*(fdr == min(fdr)), by = lab]
peaks$summit[ peaks$ClusterType == 'Singleton'] = 1

singleton_labs = peaks$lab[ peaks$ClusterType == 'Singleton']
peaks$lab[ peaks$ClusterType == 'Singleton'] = paste(singleton_labs, 1:length(singleton_labs),sep='')

peaks$bin1_end = peaks$bin1_mid + RESOLUTION
peaks$bin2_end = peaks$bin2_mid + RESOLUTION
peaks_final = subset(peaks, select = c("chr", "bin1_mid", "bin1_end", "chr", "bin2_mid", "bin2_end", "count", "expected2", "fdr", "lab", "ClusterSize", "ClusterType", "NegLog10P", "summit"))
colnames(peaks_final) = c('chr1', 'start1', 'end1', 'chr2', 'start2', 'end2', 'count', 'expected', 'fdr', 'ClusterLabel', 'ClusterSize', 'ClusterType', 'ClusterNegLog10P', 'ClusterSummit')

fout = paste(INDIR,SET,'.',FDR,'.sig3Dinteractions.bedpe',sep='')
write.table(peaks_final, fout, row.names = FALSE, col.names = TRUE, quote=FALSE, sep='\t')
