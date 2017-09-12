library(dplyr)

# Load from MotifDb, which already has the data
md.tbl <- read.table("../../MotifDb/inst/scripts/import/jaspar2016/md-raw.tsv",
                     header  = TRUE, sep = "\t")

# Make the important columns into a new df
motif.class <- data_frame(motif = paste0(md.tbl$ma.name,".",md.tbl$unknown),
                          TF = md.tbl$gene.symbol,                          
                          class = md.tbl$class,                          
                          family = md.tbl$family)


save(motif.class,file =  "./motif_class_pairs.Rdata")

                   
