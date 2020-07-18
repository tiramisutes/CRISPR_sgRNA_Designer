#!/usr/bin/Rscript

argv <- commandArgs(T)
## ------------------------------------------------------------------------
library(tidyverse)
library(iSTOP)
library("BSgenome.Gossypiumhirsutum.HZAU.Ghirsutum")
library("TxDb.GossypiumhirsutumGFF.HZAU.Ghirsutum")
genome <- BSgenome.Gossypiumhirsutum.HZAU.Ghirsutum
txdb <- TxDb.GossypiumhirsutumGFF.HZAU.Ghirsutum
mytranid <- argv[1]
## ------------------------------------------------------------------------
## Extract the cds ranges grouped by transcript from 'txdb':
cds <- cdsBy(txdb, by="tx", use.names=TRUE)
dfcds <- tibble(as.data.frame(cds)) %>% 
  dplyr::select(tx = group_name,
                gene = group_name,
                exon = exon_rank,
                chr = seqnames,
                strand = strand,
                start = start,
                end = end) %>% 
  dplyr::mutate(gene = gsub("\\..*", "", gene),
                chr = as.character(chr),
                exon = as.integer(exon),
                strand = as.character(strand))
## ------------------------------------------------------------------------
# (1) Detect iSTOP targets for a single transcripts of a gene
singlegene <- 
  dfcds %>%
  filter(tx %in% mytranid) %>%
  locate_codons(genome) %>%
  locate_PAM(genome)
#add_RFLP(width = 150) %>%  # 150 is the maximum allowed
#add_RFLP(recognizes = 't') # Enzymes that cut if c edits to t
## ------------------------------------------------------------------------
# (2) Save results to csv
write.csv(singlegene, file = argv[2], row.names = FALSE)