args = commandArgs(trailingOnly=TRUE)
infile = args[1]
outfile = args[2]
method = args[3]
fdr_thres = args[4]

tab = read.table(infile)

fdr = p.adjust(p = tab$V7, method = method)

tab$fdr = fdr 

tabout = tab[tab$fdr <= fdr_thres,]
tabout = tabout[,c(1,2,3,4,5,6,8)]

write.table(tabout,file = outfile,sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE)
