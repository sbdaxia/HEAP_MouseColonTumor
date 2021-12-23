library(tidyverse)
mir_peaks <- readRDS("Datafiles/miRNA-merged-peaks-list-09282019-withIDs.rds")

mir_194 <-  mir_peaks$`miR-194-5p` %>% as.data.frame()

write.csv(mir_194, "crc_HEAP_miR-194-5p_targets.csv", row.names = FALSE)
