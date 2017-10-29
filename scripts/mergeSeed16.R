libs <- c(
    'dplyr',
    'tidyr',
    'stringr',
    'ggplot2',
    'GenomicRanges',
    'RPostgreSQL',
    'doMC',
    'numbers',
    'doParallel'
)

for (lib in libs) {
    if (!require(lib, character.only = TRUE, quietly = TRUE)) {
        install.packages(lib, repos='http://cran.us.r-project.org')
    }
}

library(BSgenome.Hsapiens.UCSC.hg38)
hg38 <- BSgenome.Hsapiens.UCSC.hg38

source("../my_R_functions/utility_functions.R")
source("../my_R_functions/stat_functions.R")
source("../my_R_functions/plot_functions.R")
#source("/scratch/github/BDDS/trenadb/src/utils.R")
source("/scratch/github/BDDS/footprints/testdb/src/dbFunctions.R")

# Bring in the huge DF
load("/scratch/data/all.TF.fimo.samples.ratio.ALL.df.RData")

# Read data from lymphoblast
db_lymph_hint <- src_postgres(drv=dbDriver("PostgreSQL"),
                              user="trena",
                              password="trena",
                              dbname = "lymphoblast_hint_16", #Use Seed 16
                              host = "localhost")
db_lymph_well <- src_postgres(drv=dbDriver("PostgreSQL"),
                              user="trena",
                              password="trena",
                              dbname = "lymphoblast_wellington_16", # Use Seed 16
                              host = "localhost")
hint_regions <- tbl(db_lymph_hint, "regions")
hint_hits    <- tbl(db_lymph_hint, "hits")
well_regions <- tbl(db_lymph_well, "regions")
well_hits    <- tbl(db_lymph_well, "hits")

# Function for the annotation
merge_fimo_hint_wellington_one_chrom <- function(chrom_str,
                                                 fimo_tbl,
                                                 hint_regions_tbl,
                                                 hint_hits_tbl,
                                                 well_regions_tbl,
                                                 well_hits_tbl
                                                 ) {

    # some tables use chr22 and some just use 22
    chrom_long_str = paste("chr",chrom_str, sep="")

    # select one chromosome from my data table
    fimo_tbl %>%
        filter(chrom==chrom_str) %>%
        select(-empty) ->
        chrom_all_tf_df

    # select one chromosome from hint
    hint_regions_tbl %>%
        filter(chrom==chrom_long_str) %>%
        left_join(hint_hits_tbl, by="loc") %>%
        as.data.frame %>%
        as.tbl %>%
        select(start, endpos, strand, name, score1) %>%
        rename("score1"="h_score") ->
        chrom_hint_all_tbl

    # select one chromosome from wellington
    well_regions_tbl %>%
        filter(chrom==chrom_long_str) %>%
        left_join(well_hits_tbl, by="loc") %>%
        as.data.frame %>%
        as.tbl %>%
        select(start, endpos, strand, name, score1) %>%
        rename("score1"="w_score") ->
        chrom_well_all_tbl

    # keep only max hint score but count total nontrivial scores
    chrom_hint_all_tbl %>%
        group_by(start, endpos, name, strand) %>%
        mutate(h_count = n()) %>%
        group_by(start, endpos, name, strand) %>%
        mutate(h_max_score = max(h_score)) %>%
        distinct(start, endpos, name, strand, .keep_all = TRUE) %>%
        select(-h_score) ->
        chrom_hint_unique_tbl

    # keep only min wellington score but count total nontrivial scores
    chrom_well_all_tbl %>%
        group_by(start, endpos, name, strand) %>%
        mutate(w_count = n()) %>%
        group_by(start, endpos, name, strand) %>%
        mutate(w_min_score = min(w_score)) %>%
        distinct(start, endpos, name, strand, .keep_all = TRUE) %>%
        select(-w_score) ->
        chrom_well_unique_tbl
    
    # merge hint and wellington into my table
    chrom_all_tf_df %>%
        left_join(chrom_hint_unique_tbl, by=c("start", "endpos", "strand", "motifname"="name")) %>%
        left_join(chrom_well_unique_tbl, by=c("start", "endpos", "strand", "motifname"="name")) %>%
        replace_na(list(h_count=0, w_count=0, h_max_score=0, w_min_score=0)) ->
        chrom_all_tf_df_merged
    
    return(chrom_all_tf_df_merged)

}

# Run the big function using multicore process
# I added X and Y here too
chromosomes <- as.character(c(1:22, "X","Y"))

library(BiocParallel)
register(MulticoreParam(workers = 25, stop.on.error = FALSE, log = TRUE), default = TRUE)

all.TF.df.fimo.hint.well <- bplapply(chromosomes, merge_fimo_hint_wellington_one_chrom,
                   fimo_tbl = ,
                   hint_regions_tbl = ,
                   hint_hits_tbl = ,
                   well_regions_tbl = ,
                   well_hits_tbl = )

# Combine them
all.TF.df.fimo.hint.well <- bind_rows(all.TF.df.fimo.hint.well)

# Look at it
str(all.TF.df.fimo.hint.well)

# Save it
save(all.TF.df.fimo.hint.well, file = "/scratch/data/all.TF.df.fimo.hint.well.seed16.ALL.Rdata")

# Exploration can wait for now....add that later
