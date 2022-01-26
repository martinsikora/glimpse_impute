
###########################################
## generate GLIMPSE summary master table ##
###########################################

## --------------------------------------------------------------------------------
## libraries

suppressMessages(library(tidyverse))
suppressMessages(library(yaml))


## --------------------------------------------------------------------------------
## command line arguments

args <- commandArgs(trailingOnly = TRUE)
cfiles <- args[1]  ## list of config files to process
outpx <- args[2]


## --------------------------------------------------------------------------------
## parse individual summary tables

cf <- scan(cfiles, what = character())

d <- map_dfr(cf, ~{
    r1 <- read_yaml(.x)
    p <- paste(r1$out_dir, "/", r1$prefix , ".glimpse.summary.tsv", sep = "")
    r2 <- read_tsv(p) %>%
        mutate(prefix = .x, path = r1$out_dir) %>%
        select(prefix, everything())
    r2
})


## --------------------------------------------------------------------------------
## write output table and plot summaries

write_tsv(d, file = paste(outpx, ".glimpse.summary.tsv", sep = ""))

d1 <- d %>%
    select(sampleId, prefix, set, dpAvg) %>%
    pivot_wider(names_from = set,
                values_from = dpAvg) %>%
    mutate(ratio = target / all)

br <- 10^(-5:2)

pdf(paste(outpx, ".glimpse.summary.plot.pdf", sep = ""), width = 5, height = 4)

## coverage target / all SNPs
p <- ggplot(d1, aes(x = all, y = target))
p +
    geom_abline(slope = 1,
                intercept = 0,
                linetype = "dashed",
                size = 0.25) +
    geom_point(size = 1.5,
               alpha = 0.5,
               color = "royalblue4") +
    coord_equal() +
    scale_x_log10(breaks = br,
                  labels = formatC(br, format = "fg")) +
    scale_y_log10(breaks = br,
                  labels = formatC(br, format = "fg")) +
    xlab("Average coverage all SNPs") +
    ylab("Average coverage target SNPs") +
    annotation_logticks(sides = "bl",
                        size = 0.25) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(linetype = "dotted",
                                          size = 0.25,
                                          color = "grey"))

## coverage ratio target / all SNPs
p <- ggplot(d1, aes(x = ratio))
p +
    geom_density(size = 0.25,
               alpha = 0.5,
               color = "royalblue4",
               fill = "royalblue4") +
    scale_x_log10() +
    xlab("Ratio of average coverage target / all SNPS") +
    ylab("Density") +
    annotation_logticks(sides = "b",
                        size = 0.25) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(linetype = "dotted",
                                          size = 0.25,
                                          color = "grey"))

## total number above coverage cutoffs

cv <- c(0, 0.1, 0.5, 1, 5, 10)

d2 <- map_dfr(cv, ~{
    d %>%
        filter(dpAvg >= .x) %>%
        count(set) %>%
        mutate(cv = .x)
})

p <- ggplot(d2, aes(x = factor(cv),
                    y = n))
p +
    geom_hline(yintercept = 0) +
    geom_col(aes(fill = set),
             width = 0.7,
             position = position_dodge(width = 0.9)) +
    geom_text(aes(label = prettyNum(n, big.mark = ","),
                  group = set),
              size = 2,
              hjust = 0.5,
              vjust = 0,
              position = position_dodge(width = .9)) +
    scale_fill_manual(name = "SNP set",
                      values = c("steelblue3", "burlywood3")) +
    xlab("Minimum coverage") +
    ylab("Number of samples") +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(linetype = "dotted",
                                          size = 0.25,
                                          color = "grey"))

dev.off()
