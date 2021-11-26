library(CLIPanalyze)
library(tidyverse)
library(GenomicRanges)

hva.peak <- readRDS("peakdata.HVA.rds")
hva.peak <- hva.peak$peaks
hvak.peak <- readRDS("peakdata.HVAK.rds")
hvak.peak <- hvak.peak$peaks

hva.specific.peak <- hva.peak[!hva.peak %over% hvak.peak,]
merged.peak <- c(hva.specific.peak,hvak.peak)
names(merged.peak) <- paste0("peak", seq(1:length(merged.peak)))
merged.peak$name <- names(merged.peak)

## annotate the merged.peaks with target gene name
merged.peak$'target_gene' <- NA
for (i in 1:length(merged.peak)) {
  if (!is.na(merged.peak$utr3[i]) | !is.na(merged.peak$`utr3*`[i])) {
    gene_name <- unique(c(merged.peak$utr3[i],merged.peak$`utr3*`[i]))
    gene_name <- gene_name[!is.na(gene_name)] 
    merged.peak$'target_gene'[i] <- paste(unlist(gene_name), collapse = " ")
  }
  else {
    gene_name <- unique(c(merged.peak$exon[i], merged.peak$intron[i],merged.peak$utr5[i],merged.peak$`utr5*`[i]))
    gene_name <- gene_name[!is.na(gene_name)]
    if (length(gene_name) >0) {
      merged.peak$'target_gene'[i] <- paste(unlist(gene_name), collapse = " ")
    }
  }
}

merged.peak$padj <- NULL
merged.peak$log2FC <- NULL
saveRDS(merged.peak, "peaks-merged.rds")