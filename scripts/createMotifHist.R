# Load necessary libraries
library(dplyr)
library(purrr)
library(fst)

# Load the various datasets I created
df.1M <- read.fst("/scratch/data/all.TF.fimo.samples.ratio.1M.df.fst")
df.5M <- read.fst("/scratch/data/all.TF.fimo.samples.ratio.5M.df.fst")
df.10M <- read.fst("/scratch/data/all.TF.fimo.samples.ratio.10M.df.fst")
df.25M <- read.fst("/scratch/data/all.TF.fimo.samples.ratio.25M.df.fst")

# This is the giant one that takes forever...
df.full <- read.fst("/scratch/data/all.TF.fimo.samples.ratio.25M.df.fst")

# Run the workflow to get their motif frequencies
df.1M %>% select(motifname) %>% group_by(motifname) %>% summarise(Freqs.1M = n()) -> motif.hits.1M
df.5M %>% select(motifname) %>% group_by(motifname) %>% summarise(Freqs.5M = n()) -> motif.hits.5M
df.10M %>% select(motifname) %>% group_by(motifname) %>% summarise(Freqs.10M = n()) -> motif.hits.10M
df.25M %>% select(motifname) %>% group_by(motifname) %>% summarise(Freqs.25M = n()) -> motif.hits.25M
df.full %>% select(motifname) %>% group_by(motifname) %>% summarise(Freqs.Full = n()) -> all.motif.hits

# Smush together all the dataframes using purrr::reduce
all.hits <- list(all.motif.hits,
                 motif.hits.25M,
                 motif.hits.10M,
                 motif.hits.5M,
                 motif.hits.1M)
hist.df <- reduce(all.hits, full_join) %>% arrange(Freqs.Full)

# Save the histogram data frame for plotting purposes
save(hist.df, file = "../Rdata_files/histogramDF.Rdata")
