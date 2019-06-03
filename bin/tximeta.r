#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 6) {
  stop("Usage: tximeta.r <coldata> <salmon_out> <index> <tximeta> <fasta> <gtf>", call.=FALSE)
}

path = args[2]

fns = list.files(path, pattern = "quant.sf", recursive = T, full.names = T)
names = basename(dirname(fns))
coldata = read.csv(args[1])
coldata = coldata[match(names, coldata[,1]),]
coldata = cbind(files = fns, coldata)

library(tximeta)
setTximetaBFC(".")
annotation = strsplit(args[4], "\\.")[[1]]
message("tximeta annotation parameters", annotation)

if (length(annotation) !=4 )
  stop("tximeta parameter is malformed, expected: source.organism.version.release, but got: ", args[4])

source = annotation[1]
if (!source %in% c("Ensembl", "Encode"))
  stop("First element of annotation should be Ensembl or Enconde")

organism = annotation[2]
genome = annotation[3]
release = annotation[4]

formatted_gtf = paste0(organism, ".", genome, ".", release, ".gtf")
file.symlink(args[6], formatted_gtf)
makeLinkedTxome(indexDir=args[3], source=source, organism=organism, release=release, genome=genome, fasta=args[5], gtf=formatted_gtf, write=FALSE)

se <- tximeta(coldata)
gse <- summarizeToGene(se)

saveRDS(gse, file = "gse.rds")
saveRDS(se, file = "se.rds")

# Print sessioninfo to standard out
citation("tximeta")
sessionInfo()