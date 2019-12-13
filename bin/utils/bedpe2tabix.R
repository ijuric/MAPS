library(argparse)

parser <- parser <- ArgumentParser(description='Generate tabix file from bedpe-like files')
parser$add_argument('-i', '--bedpe', type="character", help='input bedpe file')
parser$add_argument('-t', '--tabix', type="character", help='output tabix filename')
parser$add_argument('-s', '--score', type="character", help='column name for the score data, default = fdr', default = "fdr")

#parser$print_help()

args <- parser$parse_args()
bedpe <- args$bedpe
tabix_filename <- args$tabix
score_column <- args$score


#print(score_column)
#epsilon <- 0.0000001
#quit()
d <- read.csv(bedpe, sep = "\t")
colnames(d)[1:6] <- c("chr1", "start1", "end1", "chr2", "start2", "end2")
#head(d)
#print(d$score)
d$score <- -log(as.numeric(as.character((d[,score_column]))))
#head(d)
#exit()
#d <- d[complete.cases(d),]
#print(d$score)
d$lr <- paste0(d$chr2, ":", d$start2, "-", d$end2, ",", d$score)
d$rl <- paste0(d$chr1, ":", d$start1, "-", d$end1, ",", d$score)
tabixl <- d[,c("chr1", "start1", "end1", "lr")]
tabixr <- d[,c("chr2", "start2", "end2", "rl")]
colnames(tabixr) <- colnames(tabixl)
tabixl$id <- seq(1, dim(tabixl)[1] * 2, 2)
tabixr$id <- seq(2, dim(tabixr)[1] * 2, 2)
tabix <- rbind(tabixl, tabixr)
tabix <- tabix[order(tabix$id),]
tabix$dir <- "."
options(scipen = 999)
print(colnames(tabix))
tabix <- tabix[order(tabix$chr1, tabix$start1),]
write.table(tabix, tabix_filename, quote = F, row.names = F, col.names = F, sep = "\t")
system(paste0("bgzip -c ", tabix_filename, " > ", tabix_filename, ".gz"))
system(paste0("tabix -p bed ", tabix_filename, ".gz"))
#quit()
#write.table(tabix, tabix_filename, quote = F, row.names = F, col.names = F, sep = "\t")
