
##################################
## summary stats for imputation ##
##################################


## --------------------------------------------------------------------------------
## libraries

suppressMessages(library(tidyverse))


## --------------------------------------------------------------------------------
## command line arguments

args <- commandArgs(trailingOnly = TRUE)
infiles <- args[-length(args)]
outfile <- args[length(args)]


## --------------------------------------------------------------------------------
## read and summarize infiles

d <- map_dfr(infiles, read_tsv)

r <- d %>%
    group_by(sampleId, set) %>%
    summarise(dpAvg = weighted.mean(dpAvg, nSnps),
              gpAvg = weighted.mean(gpAvg, nSnps),
              nSnps = sum(nSnps),
              .groups = "drop")


## --------------------------------------------------------------------------------
## writre output table

write_tsv(r, file = outfile)
