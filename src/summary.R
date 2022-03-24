
##################################
## summary stats for imputation ##
##################################


## --------------------------------------------------------------------------------
## libraries

suppressMessages(library(tidyverse))


## --------------------------------------------------------------------------------
## command line arguments

args <- commandArgs(trailingOnly = TRUE)
vcf <- args[1]
bed <- args[2]
outfile <- args[3]


## --------------------------------------------------------------------------------
## target sites

p <- pipe(paste("gzip -fcd ", bed, " | awk '{print $1\":\"$3}'", sep = ""))
targets <- scan(p, what = character())
close(p)


## --------------------------------------------------------------------------------
## annotations

## sample IDs
p <- pipe(paste("bcftools query -l ", vcf, sep = ""))
samples <- scan(p, what = character())
close(p)
nSamples <- length(samples)

## SNP IDs
p <- pipe(paste("bcftools query -f '%CHROM:%POS\n' ", vcf, sep = ""))
snps <- scan(p, what = character())
close(p)

isTarget <- which(snps %in% targets)
chrom <- snps[1] %>%
    strsplit(":") %>%
    map_chr(1)


## --------------------------------------------------------------------------------
## DP/GP average for each sample chunk

## chunk samples into sets of 50 to avoid issues with too large matrices
chunkSize <- 50
nChunks <- (nSamples - 1) %/% chunkSize + 1

samplesChunk <- split(samples, rep(1:nChunks, each = chunkSize, length.out = nSamples))

d <- map_dfr(1:nChunks, ~{

    cat("__ processing sample chunk", .x , "/", nChunks, "__\n")
    samplesStr <- paste(samplesChunk[[.x]], collapse = ",")

    ## DP average
    p <- pipe(paste("bcftools query -f '[%DP\t]\n' -s ", samplesStr, " ", vcf, sep = ""))
    r <- read_tsv(p, col_names = samplesChunk[[.x]], col_types = str_dup("d", length(samplesChunk[[.x]]) + 1)) %>%
        select(all_of(samplesChunk[[.x]])) %>%
        as.matrix() %>%
        replace_na(0)

    dpA <- colMeans(r)
    dpT <- colMeans(r[isTarget,, drop = FALSE])


    ## GP average
    p <- pipe(paste("bcftools query -f '[%GP\t]\n' -s ", samplesStr, " ", vcf, sep = ""))
    r <- read_tsv(p, col_names = samplesChunk[[.x]], col_types = str_dup("c", length(samplesChunk[[.x]]) + 1)) %>%
        select(all_of(samplesChunk[[.x]])) %>%
        as.matrix() %>%
        strsplit(",") %>%
        map_dbl(~ max(as.numeric(.x))) %>%
        matrix(ncol = length(samplesChunk[[.x]]))

    gpA <- colMeans(r)
    gpT <- colMeans(r[isTarget,, drop = FALSE])


    ## results table
    r1 <- tibble(sampleId = samplesChunk[[.x]],
                 chromosome = chrom,
                 set = "all",
                 nSnps = length(snps),
                 dpAvg = dpA,
                 gpAvg = gpA)

    r2 <- tibble(sampleId = samplesChunk[[.x]],
                 chromosome = chrom,
                 set = "target",
                 nSnps = length(isTarget),
                 dpAvg = dpT,
                 gpAvg = gpT)

    ## return dfr
    bind_rows(r1, r2)
})

d <- d %>%
    arrange(set, sampleId)

write_tsv(d, file = outfile)
