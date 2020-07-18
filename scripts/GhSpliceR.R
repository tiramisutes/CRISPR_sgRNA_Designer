#!/usr/bin/Rscript

argv <- commandArgs(T)
source("./scripts/GhSpliceR_global.R")
## ------------------------------------------------------------------------
# Example data
input_id <- argv[1]
input_pam <- argv[2]
input_enzyme = argv[3]
input_splice_site = argv[4]

# Render table with sgRNAs
GhSpliceR <- 
  generate_guides(ensembl_id = input_id, pam = input_pam)[[1]] %>%
  filter(enzyme == input_enzyme) %>%
  mutate(splice_site = paste0("splice-", splice_site)) %>%
  filter(mapply(FUN = grepl, pattern = splice_site, x = input_splice_site)) %>%
  dplyr::select(exon, splice_site, guide, pam, be_efficiency, be_dinucleotide, be_position, abe_efficiency, abe_dinucleotide, abe_position, full_motif, gene_id, id) %>%
  dplyr::rename(Exon = 1, `Splice-site` = 2, Protospacer = 3, PAM = 4,
                `CBE efficiency` = 5, `CBE dinucleotide` = 6, `CBE position` = 7,
                `ABE efficiency` = 8, `ABE dinucleotide` = 9, `ABE position` = 10,
                `Motif` = 11, `Gene ID` = 12, `Transcript ID` = 13)
write.csv(GhSpliceR, file = argv[5], row.names = FALSE)