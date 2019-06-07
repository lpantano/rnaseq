#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 6) {
  stop("Usage: tximeta.r <coldata> <salmon_out> <index> <tximeta> <fasta> <gtf>", call.=FALSE)
}

path = args[2]
coldata = args[1]

tx2gene = "tx2gene.csv"
info = file.info(tx2gene)
if (info$size == 0){
  tx2gene = NULL
}else{
  rowdata = read.csv(tx2gene, header = FALSE)
  colnames(rowdata) = c("tx", "gene_id", "gene_name")
  tx2gene = rowdata[,1:2]
}

if (!file.exists(args[5]))
    stop("FASTA file doesn't exist ", args[6])
if (!file.exists(args[6]))
    stop("GTF file doesn't exist ", args[7])

fns = list.files(path, pattern = "quant.sf", recursive = T, full.names = T)
names = basename(dirname(fns))
names(fns) = names
if (file.exists(coldata)){
    coldata = read.csv(args[1])
    coldata = coldata[match(names, coldata[,1]),]
    coldata = cbind(files = fns, coldata)
}else{
    message("ColData not avaliable ", coldata)
    coldata = data.frame(files = fns, names = names)
}

if (args[4] != "NULL"){
  annotation = strsplit(args[4], "\\.")[[1]]
  message("tximeta annotation parameters", annotation)

  if (length(annotation) !=4 )
    stop("tximeta parameter is malformed, expected: source.organism.version.release, but got: ", args[4])

  source = annotation[1]
  if (!source %in% c("Ensembl", "Encode"))
    stop("First element of annotation should be Ensembl or Enconde")
}

library(SummarizedExperiment)

if (args[4] == "NULL"){ # if not genome version is giving
  library(tximport)

  txi = tximport(fns, type = "salmon", txOut = TRUE)
  rownames(coldata) = coldata[["names"]]
  rowdata = rowdata[match(rownames(txi[[1]]), rowdata[["tx"]]),]
  se = SummarizedExperiment(assays = list(counts = txi[["counts"]],
                                          abundance = txi[["abundance"]],
                                          length = txi[["length"]]),
                            colData = DataFrame(coldata),
                            rowData = rowdata)
  if (!is.null(tx2gene)){
    gi = summarizeToGene(txi, tx2gene = tx2gene)
    growdata = unique(rowdata[,2:3])
    growdata = growdata[match(rownames(gi[[1]]), growdata[["gene_id"]]),]
    gse = SummarizedExperiment(assays = list(counts = gi[["counts"]],
                                            abundance = gi[["abundance"]],
                                            length = gi[["length"]]),
                              colData = DataFrame(coldata),
                              rowData = growdata)
  }
  
}else{ #when geneome version is giving like Ensembl.Hsapiens.GRCh38.69
  library(tximeta)
  setTximetaBFC(".")
  
  organism = annotation[2]
  genome = annotation[3]
  release = annotation[4]
  
  formatted_gtf = paste0(organism, ".", genome, ".", release, ".gtf")
  file.symlink(args[6], formatted_gtf)
  makeLinkedTxome(indexDir=args[3], source=source, organism=organism, release=release, genome=genome, fasta=args[5], gtf=formatted_gtf, write=FALSE)
  se <- tximeta(coldata)
  gse <- summarizeToGene(se)
}

if(exists("gse")){
  saveRDS(gse, file = "gse.rds")
  write.csv(assays(se)[["abundance"]], "merged_salmon_gene_tpm.csv")
  write.csv(assays(se)[["counts"]], "merged_salmon_gene_reads.csv")
}

saveRDS(se, file = "se.rds")
write.csv(assays(se)[["abundance"]], "merged_salmon_tx_tpm.csv")
write.csv(assays(se)[["counts"]], "merged_salmon_tx_reads.csv")

# Print sessioninfo to standard out
citation("tximeta")
sessionInfo()