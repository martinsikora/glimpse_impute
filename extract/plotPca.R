
###########################################
## generate GLIMPSE summary master table ##
###########################################

## --------------------------------------------------------------------------------
## libraries

suppressMessages(library(tidyverse))
suppressMessages(library(ggrepel))


## --------------------------------------------------------------------------------
## command line arguments

args <- commandArgs(trailingOnly = TRUE)
infile <- args[1]  ## list of prefixes to process
outfile <- args[2]


## --------------------------------------------------------------------------------
## read and plot

d <- read_tsv(infile)

lab <- paste("PC", 1:20, sep = "") %>%
    matrix(nrow = 2)

pdf(outfile, width = 7, height = 6)
for(i in 1:ncol(lab)){

    ## select PC subset
    d1 <- d %>%
        select(IID, lab[1, i], lab[2, i],)


    ## identify outliers for labelling
    d2 <- d1 %>%
        select(-IID) %>%
        as.matrix()

    md <- mahalanobis(d2, colMeans(d2), cov(d2))

    d3 <- d1 %>%
        mutate(md = md) %>%
        filter(min_rank(desc(md)) < 20)

    ## plot
    p <- ggplot(d1, aes_string(x = lab[1, i],
                               y = lab[2, i]))
    print(p +
        geom_point(size = 1.5,
                   alpha = 0.5,
                   color = "royalblue4") +
        geom_text_repel(aes(label = IID),
                        size = 2,
                        segment.color = "grey",
                        segment.size = 0.25,
                        data = d3,
                        max.overlaps = Inf) +
        theme_bw() +
        theme(panel.grid.minor = element_blank(),
              panel.grid.major = element_line(linetype = "dotted",
                                              size = 0.25,
                                              color = "grey")))
}
dev.off()
