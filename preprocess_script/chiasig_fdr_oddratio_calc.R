args = commandArgs(trailingOnly=TRUE)
infile = args[1]
outfile = args[2]
log_ratio = args[3]
method = args[4]
fdr_thres = args[5]

tab = read.table(infile)

fdr = p.adjust(p = tab$V7, method = method)
log_oddratio = log(tab$V10) - log(tab$V11)

tab$fdr = fdr
tab$log_oddratio = log_oddratio

tabout = tab[tab$fdr <= fdr_thres & tab$log_oddratio >= log_ratio,]

tabout = tabout[,c(1,2,3,4,5,6,12)]

write.table(tabout,file = outfile,sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE)