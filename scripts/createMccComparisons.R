# Using the 10 million data points, collect the data for plotting a range of MCC values
# over the corresponding ranges of Wellington and HINT cutoffs
library(FPML)
library(fst)
library(ggpubr)

# Load the dataset
seed.16 <- read.fst("/scratch/data/all.TF.df.fimo.hint.well.seed16.annotated.10M.fst")
seed.20 <- read.fst("/scratch/data/all.TF.df.fimo.hint.well.seed20.annotated.10M.fst")

# Glue them together to make the joined dataset
joined.df <- joinModelData(seed.16, seed.20)

# Create a function that:
# 1) Takes a HINT cutoff value that we'll vary
# 2) Runs "prepModelData" to filter according to the cutoff value
# 3) Creates the boosted model
# 4) Creates the linear models
# 5) Pulls the maxMCC values for the boosted and linear models and returns them

findMCCsForHINTCutoff <- function(HINT.cutoff){

    # Prep the data as specified by the cutoff
    prepped.data <- prepModelData(joined.df, "both", hintCutoff = HINT.cutoff)

    # Create the boosted model
    boosted.results <- buildBoostedModel(prepped.data)

    # Create the linear models
    linear.results <- buildLinearModels(prepped.data)

    # Pull out the MCC max value vector
    mcc.vector <- extractMaxMCC(boosted.results[[1]], linear.results)

    # Return just the boosted and full linear values
    return(mcc.vector[c("gradient boosted model",
                      "linear model (all regressors)")])    
}

# Create a similar function, but make it do the variation for Wellington
findMCCsForWellingtonCutoff <- function(Well.cutoff){

    # Prep the data as specified by the cutoff
    prepped.data <- prepModelData(joined.df, "both", wellCutoff = Well.cutoff)

    # Create the boosted model
    boosted.results <- buildBoostedModel(prepped.data)

    # Create the linear models
    linear.results <- buildLinearModels(prepped.data)

    # Pull out the MCC max value vector
    mcc.vector <- extractMaxMCC(boosted.results[[1]], linear.results)

    # Return just the boosted and full linear values
    return(mcc.vector[c("gradient boosted model",
                      "linear model (all regressors)")])
}

# Run HINT Cutoff script over a range of 20 values, collect into a DF, and plot
hint.values <- seq(0,10,0.5)
all.hint.list <- lapply(hint.values, findMCCsForHINTCutoff)
all.boosted <- sapply(all.hint.list, function(x) x[1])
all.linear <- sapply(all.hint.list, function(x) x[2])

hint.df <- data_frame(Threshold = rep(hint.values,2),
                      MCC = c(all.boosted,all.linear),
                      Model = c(rep("Gradient Boosted",length(hint.values)),
                                rep("Logistic Regression", length(hint.values))
                                )
                      )

hint.plot <- ggscatter(hint.df, x = "Threshold", y = "MCC", color = "Model")
                      
well.values <- seq(-5, 0, 0.25)
all.well.list <- lapply(well.values, findMCCsForWellingtonCutoff)
all.boosted <- sapply(all.well.list, function(x) x[1])
all.linear <- sapply(all.well.list, function(x) x[2])

well.df <- data_frame(Threshold = rep(well.values,2),
                      MCC = c(all.boosted,all.linear),
                      Model = c(rep("Gradient Boosted",length(well.values)),
                                rep("Logistic Regression", length(well.values))
                                )
                      )

well.plot <- ggscatter(well.df, x = "Threshold", y = "MCC", color = "Model")

# Make the figure 1 x 2

ggarrange(hint.plot, well.plot,
          ncol = 2, nrow = 1,
          common.legend= TRUE) %>%
    ggexport(filename = "MccComparisonPlots.png")

# Save the data
save(hint.df, well.df, hint.plot, well.plot, file = "./MccComparisonResults.Rdata")
