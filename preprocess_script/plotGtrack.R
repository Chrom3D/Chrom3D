args = commandArgs(trailingOnly=TRUE)
file = args[1]
selChr = args[2]
outFile = args[3]

stopifnot(grep(".png", outFile) == 1)

a = read.table(file, colClasses=c("character", "numeric","numeric", "character", "numeric", "numeric", "character" ))



sel = a[a$V1==selChr,]

minPos = min(sel$V2)
maxPos = max(sel$V3)

png(outFile, width=1800, height=1800, res=300)
#par(xpd=NA)
plot(1, type="n", xlab="", ylab="", xlim=c(minPos, maxPos), ylim=c(minPos, maxPos), xaxt="n", yaxt="n")


mi=1.02
ma=1.04

for(i in 1:nrow(sel)) {
  if(a$V7[i] != ".") {
    links = strsplit(sel$V7[i], ";")
    for(l in links[[1]]) {
      if(strsplit(l,":")[[1]][1] == selChr) {
        linkleft = as.numeric(strsplit(strsplit(l,":")[[1]][2], "-")[[1]][1])
        linkright = as.numeric(strsplit(strsplit(l,":")[[1]][[2]], "-")[[1]][2])
        rect(sel$V2[i], linkleft,  sel$V3[i], linkright, col="red", border=NA)
        rect(linkleft, sel$V2[i], linkright,sel$V3[i],  col="red",border=NA)
      }     
    }
  }  
}
par(xpd=NA)
for(i in 1:nrow(sel)) {
  if(sel$V6[i]==1) {
    rect(sel$V2[i], maxPos*mi, sel$V3[i], maxPos*ma, col="blue", border=NA)
    rect(maxPos*mi, sel$V2[i], maxPos*ma, sel$V3[i], col="blue", border=NA)
  }
}

ax=par('xaxp')
axi = seq(ax[1], ax[2], length.out=ax[3]+1)

axis(1, at=axi, labels=axi / 1e+06, las=1, cex.axis=2)
axis(2, at=axi, labels=axi / 1e+06, las=1, cex.axis=2)
dev.off()
