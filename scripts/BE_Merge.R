log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(tidyverse)   # for read_tsv()
library(purrr)  # for map(), reduce()

files <- snakemake@input[[1]]    # get file names

data <- files %>%
  # read in all the files, appending the path before the file name
  map(~ read_csv(.)) %>% 
  reduce(rbind)

write_csv(data, snakemake@output[[1]])
