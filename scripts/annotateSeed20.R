libs <- c(
    'tidyverse',
    'stringr',
    'GenomicRanges',
    'RPostgreSQL',
    'doMC',
    'numbers',
    'doParallel',
    'Matrix',
    'fst'
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
all.TF.df.fimo.hint.well <-
    read.fst(path = "/scratch/data/all.TF.df.fimo.hint.well.seed20.10M.Rdata")

# Filter the motif class using the dataframe
filtered.motif.class <- semi_join(fixed.motif.class,
                                  all.TF.df.fimo.hint.well,
                                  by = "motifname")

# Create the one-hot motif class
filtered.motif.class %>%
    # clean up and subset to only relevant motifs
    mutate_all(str_trim) %>%
    # fix double classes
    mutate(class = str_split(class, "::")) %>%
    unnest(class) %>%
    # create one-hot(ish, some double matches) version
    mutate(dummy_yesno = 1) %>%
    distinct %>%
    spread(class, dummy_yesno, fill = 0) ->
    motif_class_hot

# Add the motif class to the dataframe
all.TF.df.fimo.hint.well.annotated <- left_join(all.TF.df.fimo.hint.well,
                                                motif_class_hot)

# Add GC content
get_gc_content <- function(start_col, end_col, chrom_col, b=100) {
    require(GenomicRanges)

    window_center <- round((start_col + end_col)/2)
    windows <- getSeq(hg38, paste0("chr",chrom_col), window_center-b, window_center+b)

    alph_freq <- alphabetFrequency(windows)
    gc_content <- rowSums(alph_freq[,c("C","G")])/(2*b+1)

    return(gc_content)
}


all.TF.df.fimo.hint.well.annotated %>%
    mutate("gc_content" = get_gc_content(start,endpos,chrom)) ->
    all.TF.df.fimo.hint.well.annotated

# Add TSS Distance
# Changed host to localhost
db_gtf <- dbConnect(PostgreSQL(), user= "trena", password="trena", dbname="hg38", host="localhost")
query <- "select * from gtf where moleculetype='gene' and gene_biotype='protein_coding'"
tss_raw_table <- dbGetQuery(db_gtf, query)[, c("chr", "gene_name", "start", "endpos","strand")]

tss_raw_table %>%
    mutate(ref = ifelse(strand == '+', start, endpos)) %>%
    select("chrom" = "chr", "ts_start" = "ref") %>%
    filter(chrom != 'chrMT') %>%
    mutate(chrom=str_sub(chrom,  start = 4)) ->
    tss_tbl

motif_gr <- makeGRangesFromDataFrame(all.TF.df.fimo.hint.well.annotated, start.field="start", end.field="endpos")
tss_gr <- makeGRangesFromDataFrame(tss_tbl, start.field="ts_start", end.field="ts_start")
dist_to_nearest_tss <- distanceToNearest(motif_gr, tss_gr, select="arbitrary")
tss_dists <- mcols(dist_to_nearest_tss)[,1]

all.TF.df.fimo.hint.well.annotated %>%
    mutate(asinh_tss_dist = asinh(tss_dists)) ->
        all.TF.df.fimo.hint.well.annotated

# Change counts to fracs
all.TF.df.fimo.hint.well.annotated %>%
    mutate(h_frac = h_count/max(h_count)) %>%
    mutate(w_frac = w_count/max(w_count)) %>%
    select(-one_of("h_count","w_count")) %>%
    select(motifname:w_min_score, h_frac, w_frac, gc_content, asinh_tss_dist, everything()) ->
    all.TF.df.fimo.hint.well.annotated

# Save the data to an fst file
write.fst(all.TF.df.fimo.hint.well.annotated,
          path="/scratch/data/all.TF.df.fimo.hint.well.seed20.annotated.10M.fst")





