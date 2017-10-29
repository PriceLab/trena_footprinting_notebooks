libs <- c(
    'tidyverse',
    'stringr',
    'GenomicRanges',
    'RPostgreSQL',
    'doMC',
    'numbers',
    'doParallel',
    'Matrix'

)

for (lib in libs) {
    if (!require(lib, character.only = TRUE, quietly = TRUE)) {
        install.packages(lib, repos='http://cran.us.r-project.org')
    }
}

library(BSgenome.Hsapiens.UCSC.hg38)
hg38 = BSgenome.Hsapiens.UCSC.hg38

source("../my_R_functions/utility_functions.R")
source("../my_R_functions/stat_functions.R")
source("../my_R_functions/plot_functions.R")
source("/scratch/github/BDDS/footprints/testdb/src/dbFunctions.R")

# Load data
load("../Rdata_files/Tfmotifmap.Rdata")
load("../Rdata_files/motif_class_pairs.Rdata")

# Make Jaspar translation table
suppressMessages(library(MotifDb))
jaspar.motifs <- subset(MotifDb, dataSource == "jaspar2016")
jaspar.df <- data_frame(Long.Name = names(jaspar.motifs),
                        Short.Name = trimws(values(jaspar.motifs)$providerName)
                        )

fixed.motif.class <- motif.class %>%
    left_join(jaspar.df, by = c("motif" = "Short.Name")) %>%
    select("motifname" = "Long.Name", class)

# Grab the FIMO dataframe (~400 million records)
load("/scratch/data/all.TF.df.fimo.hint.well.seed20.ALL.Rdata")
