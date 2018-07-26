
## run example:
## /opt/R-3.4.3/lib64/R/bin/Rscript MAPS_regression.r /home/jurici/work/PLACseq/MAPS2/results/bing_mESC_intersect_subsamples/ 
## MY_113.MY_115 19 RH_129-130.uniq.paired.sorted.nodup.nsrt.5k.MAPS2_filter
##
## arguments:
## INFDIR - dir with reg files
## SET - dataset name
## chroms - number of chromosomes (19 for mouse, 22 for human)
## FILTER - file containing bins that need to be filtered out. Format: two columns "chrom", "bin". "chrom" contains 'chr1','chr2',.. "bin" is bin label


library(VGAM)
options(warn=1)
 
args <- commandArgs(trailingOnly=TRUE)

fltr = data.frame(chrom='chrNONE',bin=-1)

if (length(args) != 4) {
    print('Wrong number of arguments. Stopping.')
    print('Arguments needed (in this order): INFDIR, SET, chroms, FILTER.')
    print('FILTER is optional argument. Omitt it if no filtering required.')
    print(paste('Number of arguments entered:',length(args)))
    print('Arguments entered:')
    print(args)
    quit()
} else {
    print(args)
    INFDIR = args[1]
    SET = args[2]
    chroms = paste('chr',seq(1,as.integer(args[3]),1),sep='')
    if (args[4] != 'None') {
        FILTER = args[4]
        fltr = read.table(FILTER,header=T)
    }
}

file_length_xor = 0
file_length_and = 0
for (i in chroms) {
    for (j in c('.and','.xor')) {
        print(paste('doing chromosome ',i,sep=''))
        inf_name = paste(INFDIR,'reg.',i,'.',SET,j,sep='')
        print(inf_name)
        mm = read.table(inf_name,header=T)
        if (length(mm$bin1_mid > 0)) {
            mm = mm [ abs(mm$bin1_mid - mm$bin2_mid) > 1,]  ## removing adjacent bins
            mm$chrom = i
            mm = subset(mm, !(mm$chrom %in% fltr$chrom & mm$bin1_mid %in% fltr$bin | mm$bin2_mid %in% fltr$bin )) ## filtering out bad bins
        }
        if (length(mm$bin1_mid > 0)) {
            if (j == '.and') {
                file_length_and = file_length_and + length(mm$bin1_mid)
            } else {
                file_length_xor = file_length_xor + length(mm$bin1_mid)
            }
        }
    }
}

file_length = file_length_and + file_length_xor
for (i in chroms) {
    for (j in c('.and','.xor')) {
        print(paste('doing chromosome ',i,sep=''))
        inf_name = paste(INFDIR,'reg.',i,'.',SET,j,sep='')
        outf_name = paste(inf_name,'.MAPS_pospoisson',sep = '')
        mm = read.table(inf_name,header=T)
        mm$chrom = i
        mm = mm [ abs(mm$bin1_mid - mm$bin2_mid) > 1,]  ## removing adjacent bins
        mm = subset(mm, !(mm$chrom %in% fltr$chrom & mm$bin1_mid %in% fltr$bin | mm$bin2_mid %in% fltr$bin )) ## filtering out bad bins
        fit <- vglm(count ~ logl + loggc + logm + logdist + logShortCount, family = pospoisson(), data = mm)
        mm$expected = fitted(fit)
        mm$p_val = ppois(mm$count, mm$expected, lower.tail = FALSE, log.p = FALSE) / ppois(0, mm$expected, lower.tail = FALSE, log.p = FALSE)
        m1 = mm[ mm$p_val > 1/length(mm$p_val),]
        fit <- vglm(count ~ logl + loggc + logm + logdist + logShortCount, family = pospoisson(), data = m1)
        coeff<-round(coef(fit),10)
        #m1$expected2 = fitted(fit)
        mm$expected2 <- round(exp(coeff[1] + coeff[2]*mm$logl + coeff[3]*mm$loggc + coeff[4]*mm$logm + coeff[5]*mm$logdist + coeff[6]*mm$logShortCount), 10)
        mm$expected2 <- mm$expected2 /(1-exp(-mm$expected2))
        mm$ratio2 <- mm$count / mm$expected2
        mm$p_val_reg2 = ppois(mm$count, mm$expected2, lower.tail = FALSE, log.p = FALSE) / ppois(0, mm$expected2, lower.tail = FALSE, log.p = FALSE)
        mm$p_bonferroni = mm$p_val_reg2 * file_length
        mm$fdr <- p.adjust(mm$p_val_reg2, method='fdr')
        write.table(mm,outf_name,row.names = TRUE,col.names = TRUE,quote=FALSE)
    }
}

