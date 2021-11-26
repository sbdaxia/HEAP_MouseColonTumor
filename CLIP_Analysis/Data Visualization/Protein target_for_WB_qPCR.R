library(readr)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(EnsDb.Mmusculus.v79)
library(tidyverse)


weird_protein <- read_csv("weird HVAK target proteins.csv")

annotations_org <- AnnotationDbi::select(org.Mm.eg.db,
                                         keys = as.character(weird_protein$Protein.Id),
                                         columns = c("ENSEMBL", "SYMBOL"),
                                         keytype = "UNIPROT")

non_duplicates_idx <- which(duplicated(annotations_org$ENSEMBL) == FALSE)

annotations_org <- annotations_org[non_duplicates_idx, ]

is.na(annotations_org$ENSEMBL) %>%
  which() %>%
  length()

annotations_org <- annotations_org[!is.na(annotations_org$ENSEMBL), ]

# load the RNA-Seq data
dif_rna <- read_csv("~/OneDrive - Harvard University/Haigis Lab/Projects/Halo-Ago2/Halo-Ago-KRas/Raw Data/RNA-Seq/Mouse colon tumor/4-OHT enema model/Analysis/Differential Analysis_filtered.csv")

overlap <- intersect(annotations_org$ENSEMBL, dif_rna$X1)

target <- dif_rna %>% filter(X1 %in% overlap) %>% filter(log2FoldChange > 0) %>% filter(padj < 0.05) %>% select(X1)

target_uniprot <- inner_join(target, annotations_org, by = c("X1" = "ENSEMBL"))

write_csv(target_uniprot, "Targets for WB.csv")
