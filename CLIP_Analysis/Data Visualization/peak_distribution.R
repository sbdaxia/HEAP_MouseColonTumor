library(data.table)
library(tidyverse)
library(rtracklayer)
library(Rsamtools)
library(DESeq2)
library(Rsubread)
library(CLIPanalyze)
library(RColorBrewer)
library(gplots)
library(Biostrings)

peak.data.hva<- readRDS("~/OneDrive - Harvard University/Haigis Lab/Projects/Halo-Ago2/Halo-Ago-KRas/Raw Data/HEAP-CLIP/HEAP_09282019/Analysis/CLIP_Analysis/Data Visualization/peakdata.HVA.rds")
peaks.all.hva <- peak.data.hva$peaks
peaks.all.hva <- subset(peaks.all.hva, width > 20)
peaks.all.hva <- subset(peaks.all.hva, log2FC > 0)

count.threshold <- 10
norm.counts <- counts(peak.data.hva$peak.counts, normalized = TRUE)
norm.counts <- norm.counts[names(peaks.all.hva), ]
colnames(norm.counts)[1:12] <- c(paste0("HVA", 1:6), paste0("HVA-IC", 1:6))
selected.peaks <- rowMeans(norm.counts[, paste0("HVA", 1:6)]) > count.threshold 
peaks.all.hva <- peaks.all.hva[selected.peaks, ]

padj.threshold <- 5*1e-2
hva.peaks <- subset(peaks.all.hva, padj < padj.threshold)
length(hva.peaks)
table(hva.peaks$annot)

peak.data.hvak<- readRDS("~/OneDrive - Harvard University/Haigis Lab/Projects/Halo-Ago2/Halo-Ago-KRas/Raw Data/HEAP-CLIP/HEAP_09282019/Analysis/CLIP_Analysis/Data Visualization/peakdata.HVAK.rds")
peaks.all.hvak <- peak.data.hvak$peaks
peaks.all.hvak <- subset(peaks.all.hvak, width > 20)
peaks.all.hvak <- subset(peaks.all.hvak, log2FC > 0)

norm.counts <- DESeq2::counts(peak.data.hvak$peak.counts, normalized = TRUE)
norm.counts <- norm.counts[names(peaks.all.hvak), ]
colnames(norm.counts)[1:12] <- c(paste0("HVAK", 1:6), paste0("HVAK-IC", 1:6))
selected.peaks <- rowMeans(norm.counts[, paste0("HVAK", 1:6)]) > count.threshold 
peaks.all.hvak <- peaks.all.hvak[selected.peaks, ]

hvak.peaks <- subset(peaks.all.hvak, padj < padj.threshold)
length(hvak.peaks)
table(hvak.peaks$annot)


# Make staked bargraphs for peak fractions
library(ggplot2)
library(tidyverse)
library(reshape2)
hva_peak_frac <- read_csv("feature_frac_hva_peaks.csv")
hvak_peak_frac <- read_csv("feature_frac_hvak_peaks.csv")
hva_peak_frac_df <- melt(hva_peak_frac)
colnames(hva_peak_frac_df) <- c("region","rank","fraction")
hva_peak_frac_df <- hva_peak_frac_df %>% filter(!region %in% c("utr3*","utr5*"))
hvak_peak_frac_df <- melt(hvak_peak_frac)
colnames(hvak_peak_frac_df) <- c("region","rank","fraction")
hvak_peak_frac_df <- hvak_peak_frac_df %>% filter(!region %in% c("utr3*","utr5*"))

pdf("PDF_figure/HVA_peak_genomic_fraction.pdf",
    height = 4,
    width = 6)
ggplot(hva_peak_frac_df, aes(fill=region, y=fraction, x=rank)) + 
  geom_bar(position="fill", stat="identity")+
  xlab("Top n peaks ranked by adjusted p-vaue") +
  ylab("Distribution of peaks\nacross genomic annotations") +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90),
        axis.line.x = element_line(size = 0.5),
        axis.line.y = element_line(size = 0.5)) + 
  scale_y_continuous(expand = c(0, 0))
dev.off()

pdf("PDF_figure/HVAK_peak_genomic_fraction.pdf",
    height = 4,
    width = 6)
ggplot(hvak_peak_frac_df, aes(fill=region, y=fraction, x=rank)) + 
  geom_bar(position="fill", stat="identity")+
  xlab("Top n peaks ranked by adjusted p-vaue") +
  ylab("Distribution of peaks\nacross genomic annotations") +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90),
        axis.line.x = element_line(size = 0.5),
        axis.line.y = element_line(size = 0.5)) + 
  scale_y_continuous(expand = c(0, 0))
dev.off()