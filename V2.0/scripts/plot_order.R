#usr/bin/R

args=(commandArgs(TRUE))

library("plotrix", lib.loc=args[3])
matrix<-read.table(args[1], header=TRUE, sep="\t")
numrow<-nrow(matrix)
factor<-log(numrow, numrow)
factor
ystart = 1.01

pdf(file=args[2], width=15, height=factor*100, pointsize=factor*12, onefile=TRUE)
par(mar=c(0,0,0,0))
plot.new()
boxed.labels(0.98, 1, "ATP", bg="yellow", border=TRUE, xpad=1.5, ypad=1.5, cex=1)
boxed.labels(0.98, 0.997, "COX", bg="dark green", border=TRUE, xpad=1.5, ypad=1.5, cex=1) #0.997
boxed.labels(0.98, 0.994, "COB", bg="purple", border=TRUE, xpad=1.5, ypad=1.5, cex=1) #0.994
boxed.labels(0.98, 0.991, "NAD", bg="dark blue", border=TRUE, xpad=1.5, ypad=1.5, cex=1) #0.991
boxed.labels(0.98, 0.988, "rRNA", bg="dark orange", border=TRUE, xpad=1.5, ypad=1.5, cex=1) #0.988
boxed.labels(0.98, 0.985, "tRNA", bg="dark red", border=TRUE, xpad=1.5, ypad=1.5, cex=1) #0.985
#boxed.labels(0.98, 0.982, "CR", bg="grey", border=TRUE, xpad=1.5, ypad=1.5, cex=1)

for (row in 1:nrow(matrix)) {
  pattern_num <- as.character(matrix[row,1])
  total <- matrix[row,2]
  gene_order <- as.character(matrix[row,-(1:2)])
  gene_order <- strsplit(gene_order, " ")[[1]]
  len <- length(gene_order)
  ystart <- ystart - factor*0.008 # 0.008
  linestart = ystart + factor*0.001 # 0.001
  lines(c(0,0.95), c(linestart,linestart), col="dark grey", lty=2)
  string <- c(pattern_num, "(",total, "species )")
  textbox(c(0,1), ystart, string, justify='l', margin=-0.005, box=FALSE, cex=1) # cex=1
  xstart = 0.02
  ystartfwd = ystart - factor*0.003 # 0.003
  ystartrev = ystartfwd - factor*0.002 # 0.002
  for (col in 1:len) {
    gene <- gene_order[col]
    if (grepl("COX", gene)) {
      color <- c("dark green")
      gene <- gsub("COX", "", gene)
    } else if (grepl("NAD", gene)) {
      color <- c("dark blue")
      gene <- gsub("NAD", "", gene)
      gene <- gsub("4L", "4l", gene)
    } else if (grepl("ATP", gene)) {
      color <- c("yellow")
      gene <- gsub("ATP", "", gene)
    } else if (grepl("rrn", gene)) {
      color <- c("dark orange")
      gene <- gsub("rrn", "", gene)
    } else if (grepl("trn", gene)) {
      color <- c("dark red")
      gene <- gsub("trn", "", gene)
    } else if (grepl("COB", gene)) {
      color <- c("purple")
      gene <- gsub("COB", "C", gene)
    } else {
      color <- c("grey")
      gene <- gsub("CR", "cr", gene)
    }
    if (!grepl("^-", gene)) {
       boxed.labels(xstart, ystartfwd, gene, bg=color, border=TRUE, xpad=1.5, ypad=1.5, cex=1) # cex=1
    } else {
      gene <- gsub("-", "", gene)
      boxed.labels(xstart, ystartrev, gene, bg=color, border=TRUE, xpad=1.5, ypad=1.5, cex=1) # cex=1
    }
    xstart <- xstart + factor*0.015 #0.015
  }
}
dev.off()
