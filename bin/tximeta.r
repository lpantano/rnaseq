#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 9) {
  stop("Usage: tximeta.r <salmon_out> <index> <source> <organism> <release> <genome> <fasta> <gtf>", call.=FALSE)
}

path = args[2]

fns = list.files(path, pattern = "quant.sf", recursive = T, full.names = T)
names = basename(dirname(fns))
coldata = read.csv(args[1])
coldata = coldata[match(names, coldata[,1]),]
coldata = cbind(files = fns, coldata)

library(tximeta)
setTximetaBFC(".")
formatted_gtf = paste0(args[5], ".", args[7], ".", args[6], ".gtf")
file.symlink(args[9], formatted_gtf)
makeLinkedTxome(indexDir=args[3], source=args[4], organism=args[5], release=args[6], genome=args[7], fasta=args[8], gtf=formatted_gtf, write=FALSE)


se <- tximeta(coldata)
gse <- summarizeToGene(se)

saveRDS(gse, file = "gse.rds")
saveRDS(se, file = "se.rds")

# Print sessioninfo to standard out
citation("tximeta")
sessionInfo()