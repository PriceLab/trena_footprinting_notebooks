library(RPostgreSQL)
library(dplyr)
library(GenomicRanges)
library(doParallel)
library(dbplyr)
library(BiocParallel)


source("../../BDDS/footprints/testdb/src/dbFunctions.R")

# Load the chipseq data locally
load("../Rdata_files/chipSeqData.Rdata")

# Load the TF-Motif mapping data
load("../Rdata_files/Tfmotifmap.Rdata")

# Grab all motifs for the TFs we have
allmots <- c()

for (TFname in names(TFs.to.motifs)) {
    allmots  <-  c(allmots, TFs.to.motifs[[TFname]])
}
length(unique(allmots))

# Define the function to pull motifs for a TF
create.TF.df <- function(TF, verbose = FALSE){

    # Make the database connection
    db.fimo.dplyr <- src_postgres(drv = dbDriver("PostgreSQL"),
                                  user = "trena",
                                  password = "trena",
                                  dbname = "fimo",
                                  host = "localhost")
    tbl.fimo.dplyr <- tbl(db.fimo.dplyr, "fimo_hg38")

    # Grab all hits for a TF, then grab the regions of those hits
    chipseq.hits.TF <- subset(chipseq.hits, name == TF)
    locs.TF <- chipseq.hits.TF$loc
    chipseq.regions.TF <- subset(chipseq.regions, loc %in% locs.TF)

    # this is the slow step -- doing SQL queries on tbl.fimo.dplyr = call to whole fimo database
    # need branch since %in% conversion to SQL doesn't work on length == 1
    ## Basically: we find all instances of the TF's motif(s) in FIMO
    if (length(TFs.to.motifs[[TF]]) > 1 ) {
        fimo.motifs.for.TF <- as.tbl(as.data.frame(filter(tbl.fimo.dplyr, motifname %in% TFs.to.motifs[[TF]])))
    } else {
        fimo.motifs.for.TF <- as.tbl(as.data.frame(filter(tbl.fimo.dplyr, motifname  ==  TFs.to.motifs[[TF]])))          
    }
    # Print if requested
    if (verbose == TRUE) {
        if (length(TFs.to.motifs[[TF]])==1) {
            message(paste(TF, "- querying fimo database for", length(TFs.to.motifs[[TF]]), "motif"))
        } else {
            message(paste(TF, "- querying fimo database for", length(TFs.to.motifs[[TF]]), "motifs"))
        }
    }

    # find intersect using fast genomic ranges data structure
    # We make GR objects for the fimo motifs we just found and for the chipseq regions, then find their overlaps
    gr.fimo.TF <- with(fimo.motifs.for.TF, GRanges(chrom, IRanges(start=start, end=endpos)))
    gr.chipseq.TF <- with(chipseq.regions.TF, GRanges(chrom, IRanges(start=start, end=endpos)))
    overlaps.gr.TF <- findOverlaps(gr.chipseq.TF, gr.fimo.TF, type="any")
    overlaps.TF <- as.tbl(as.data.frame(overlaps.gr.TF))

    # row numbers in fimo.motifs.for.TF where motifs overlap with chipseq peaks
    positive.fimo.examples.rows.TF <- unique(overlaps.TF$subjectHits)
    positive.fimo.examples.TF.df <- fimo.motifs.for.TF[positive.fimo.examples.rows.TF,]

    # Simply take the other rows as the negative fimo examples
    negative.fimo.examples.TF.df <- dplyr::setdiff(fimo.motifs.for.TF, positive.fimo.examples.TF.df)    
        
    # annotate and collect all samples
    positive.fimo.examples.TF.df <- as.tbl(cbind(positive.fimo.examples.TF.df, "cs_hit"=1))
    negative.fimo.examples.TF.df <- as.tbl(cbind(negative.fimo.examples.TF.df, "cs_hit"=0))
    all.fimo.examples.TF.df <- as.tbl(rbind(positive.fimo.examples.TF.df,negative.fimo.examples.TF.df))

    return(all.fimo.examples.TF.df)

}

# Run it in parallel
sorted.TF.names <- sort(names(TFs.to.motifs))
N.TF <- length(sorted.TF.names)

register(MulticoreParam(workers = 30, stop.on.error = FALSE, log = TRUE), default = TRUE)
system.time(all.TF.df <- bplapply(sorted.TF.names, create.TF.df, verbose = TRUE))
all.TF.df <- bind_rows(all.TF.df)

# Save it
fname=paste0("/scratch/data/all.TF.fimo.samples.ratio.ALL.df.RData")
save(all.TF.df, file=fname)

# Inspect it
str(all.TF.df)
