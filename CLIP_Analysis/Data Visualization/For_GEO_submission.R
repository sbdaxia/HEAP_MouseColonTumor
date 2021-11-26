library(tidyverse)

peak_info <- readRDS("Datafiles/merged-peaks-mirs-200-09282019-withID.rds")
write_csv(as.data.frame(peak_info), "Peak_information.csv")

load("../Merged_Analysis/merged_peak_analysis.rda")
write.csv(as.data.frame(hvak.hva.res), "Peak_differential_analysis.csv", row.names = TRUE)
write.csv(as.data.frame(DESeq2::counts(dds.peaks.hva.hvak, normalize = TRUE)), "Peak_normalized_counts.csv", row.names = TRUE)


