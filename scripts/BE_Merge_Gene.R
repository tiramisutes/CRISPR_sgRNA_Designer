log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(tidyverse)   # for read_tsv()
library(purrr)  # for map(), reduce()

data_path <- snakemake@input[[1]]                           # path to the data
files <- dir(data_path, pattern = snakemake@params[[1]])    # get file names

data <- files %>%
  # read in all the files, appending the path before the file name
  map(~ read_tsv(file.path(data_path, .), quote=" ")) %>% 
  reduce(rbind)

write_csv(data, snakemake@output[[1]])
