libs <- c(
    'RColorBrewer',
    'ggplot2',
    'xgboost',
    'glmnet',
    'dplyr',
    'tidyr',
    'pROC',
    'ROCR',
    'stringr',
    'caret',
    'caTools',
    'fst'
)

for (lib in libs) {
    if (!require(lib, character.only = TRUE, quietly = TRUE)) {
        install.packages(lib, repos='http://cran.us.r-project.org')
    }
}

source("../my_R_functions/utility_functions.R")
source("../my_R_functions/stat_functions.R")
source("../my_R_functions/plot_functions.R")
source("/scratch/github/BDDS/footprints/testdb/src/dbFunctions.R")

# Load motifs
load("../Rdata_files/motif_class_pairs.Rdata")

# Load the 2 datasets 
all.TF.df.fimo.hint.well.annotated <- read.fst("/scratch/data/all.TF.df.fimo.hint.well.seed20.annotated.10M.fst")

# Fix columns and filter
colnames(all.TF.df.fimo.hint.well.annotated) <- make.names(colnames(all.TF.df.fimo.hint.well.annotated), unique=TRUE)

all.TF.df.fimo.hint.well.annotated %>%
    filter(h_max_score > 3.6 |  w_min_score < -2.3) ->
    df_only_footprint_hits

# Split up the data into training/testing/validation
cols_to_drop <- c('motifname', 'chrom', 'start', 'endpos', 'strand', 'pval', 'sequence','loc')

df_only_footprint_hits %>%
    filter(chrom %in% c("2","4")) %>%
    select(-one_of(cols_to_drop)) ->
    val_df

df_only_footprint_hits %>%
    filter(chrom %in% c("1","3","5")) %>%
    select(-one_of(cols_to_drop)) ->
    test_df

df_only_footprint_hits %>%
    filter(!(chrom %in% c("1","2","3","4","5"))) %>%
    select(-one_of(cols_to_drop)) ->
    train_df

remove(all.TF.df.fimo.hint.well.annotated, df_only_footprint_hits)

# Split predictors and responses
val_df %>%
    select(-cs_hit) %>%
    as.matrix ->
    X_val

val_df %>%
    select(cs_hit) %>%
    as.matrix ->
    y_val

test_df %>%
    select(-cs_hit) %>%
    as.matrix ->
    X_test

test_df %>%
    select(cs_hit) %>%
    as.matrix ->
    y_test

train_df %>%
    select(-cs_hit) %>%
    as.matrix ->
    X_train

train_df %>%
    select(cs_hit) %>%
    as.matrix ->
    y_train

remove(val_df, test_df, train_df)

# Train the boosted model
param <- list("objective" = "binary:logistic",
              "max.depth" = 7,
              "eta" = 0.005,
              "eval.metric" = "auc"
              )

gbdt_medium <- xgboost(
    params = param,
    data = X_train,
    label = y_train,
    nround = 200,
    verbose = FALSE,
    missing = NA
)
gbdt_medium$Model.Name <- "trees with classes"

# Save the model
xgb.save(gbdt_medium,
         "../saved_models/xgboost_TF_site_predict_seed20.10M.filtered.model")

# Make the importance matrix figure
motif.class$class <- lapply(motif.class$class, make.names, unique=TRUE)
importance_matrix <- xgb.importance(colnames(X_train),model=gbdt_medium)

df <- as_data_frame(importance_matrix)
df.tf <- subset(df, Feature %in% unique(motif.class$class))
df.notf <- subset(df, !(Feature %in% unique(motif.class$class)))
tfclass.row <- c("TF_class", unname(as.list(colSums(df.tf[!(colnames(df.tf) %in% c("Feature"))]))) )
names(tfclass.row) <- colnames(df)
df.sum <- rbind(df.notf,tfclass.row)

png("../figures/seed20ImpMatrix.filtered.10M.png")
ggplot(data=df.sum, aes(x=reorder(Feature, Gain), y=Gain)) +
    geom_bar(stat="identity") +
    coord_flip() +
    theme_minimal(base_size = 30) +
    labs(x = "Feature", y="Gain")
dev.off()

# Gather stats from boosted model
medium_pred_df <- make.pred.df.from.model(gbdt_medium, X_test, y_test)
colnames(medium_pred_df)[1] <- "ChIPseq.bound"
medium_stat_df <- make.stats.df.from.preds(medium_pred_df)

# Make data ready for linear models
X_train_lin <- X_train
y_train_lin <- y_train
X_test_lin  <- X_test
y_test_lin  <- y_test

# Make the Full linear model
tf.regressors <- colnames(X_train_lin)[colnames(X_train_lin) %in% unique(motif.class$class)]
non.tf.regressors <-  colnames(X_train_lin)[!colnames(X_train_lin) %in% unique(motif.class$class)]

tf.regressors.formula <- paste("as.factor(", paste(tf.regressors, collapse=") + as.factor("), ")")
non.tf.regressors.formula <- paste(non.tf.regressors, collapse=" + ")
all.regressors.formula <- paste(non.tf.regressors.formula, tf.regressors.formula, sep=" + ")

glm.formula <- paste("ChIPseq.bound ~ ", all.regressors.formula, sep='')

glm.df.train <- as.data.frame(cbind(y_train_lin, X_train_lin)) %>% rename("cs_hit" = "ChIPseq.bound")
glm.df.test <-  as.data.frame(cbind(y_test_lin, X_test_lin)) %>% rename("cs_hit" = "ChIPseq.bound")

glm.all <- glm(as.formula(glm.formula), data=glm.df.train, family=binomial)
glm.all$Model.Name <- "glm small"

# Gather full linear stats
glm.pred.df <- make.pred.df.from.glm(glm.all, glm.df.test)
glm.stat.df <- make.stats.df.from.preds(glm.pred.df)

# Make individual linear models

stats.regressors.df <- data.frame()

for (this.regressor in colnames(X_train_lin)) {

    if (this.regressor %in% unique(motif.class$class)) {
        glm.formula <- paste("ChIPseq.bound ~ ", "as.factor(", this.regressor, ")",sep='')
    } else {
        glm.formula <- paste("ChIPseq.bound ~ ", this.regressor,sep='')
    }

    glm.df.train <- as.data.frame(cbind(y_train_lin, X_train_lin)) %>% rename("cs_hit" = "ChIPseq.bound")
    glm.df.test <-  as.data.frame(cbind(y_test_lin, X_test_lin)) %>% rename("cs_hit" = "ChIPseq.bound")

    glm.single <- glm(as.formula(glm.formula), data=glm.df.train, family=binomial)
    glm.single$Model.Name <- paste("glm ", this.regressor, sep='')

    glm.pred.single.df <- make.pred.df.from.glm(glm.single, glm.df.test)
    glm.stat.single.df <- make.stats.df.from.preds(glm.pred.single.df)

    stats.regressors.df <- rbind(stats.regressors.df, glm.stat.single.df)

}

# Gather all stats for plotting
stats.regressors.df %>%
    filter(Model.Name %in% c("glm motifscore",
                             "glm h_frac",
                             "glm h_max_score",
                             "glm w_frac",
                             "glm w_min_score",
                             "glm gc_content",
                             "glm asinh_tss_dist")) ->
    stats.regressors.filtered.df

all.stats.df <- rbind(
    medium_stat_df,
    glm.stat.df,
    stats.regressors.filtered.df)

# Save them for later
save(all.stats.df, file = "/scratch/data/seed20AllStats.10M.filtered.Rdata")

# Make a copy to play with
all.stats.orig <- all.stats.df

# Grab only the models we want for plotting
# Screen out some models for visualization
all.stats.orig %>%
    filter(Model.Name %in% c("glm h_max_score",
                             "glm w_min_score",
                             "trees with classes",
                             "glm small"
                             )) -> all.stats.df

# Change labels and print out max MCC
all.stats.df$Model.Name[all.stats.df$Model.Name== 'glm small'] <- 'linear model (all regressors)'
all.stats.df$Model.Name[all.stats.df$Model.Name== 'trees with classes'] <- 'gradient boosted model'
all.stats.df$Model.Name[all.stats.df$Model.Name== 'glm h_max_score'] <- 'HINT best score seed 20'
all.stats.df$Model.Name[all.stats.df$Model.Name== 'glm w_min_score'] <- 'Wellington best score seed 20'

max(all.stats.df$MattCC)

# Create custom color palette (matches others, mostly seed 16)
myPalette <- c("black","red","blue","green2","darkorange","purple","magenta","brown","cyan")
myColors <- myPalette[1:4]
names(myColors) <- c('gradient boosted model',
                     'linear model (all regressors)',
                     'HINT best score seed 20',
                     'Wellington best score seed 20'
                     )
colScale <- scale_colour_manual(name = "Model.Name",values = myColors)

# Reorder model names
all.stats.df$Model.Name <- factor(all.stats.df$Model.Name,
                                  levels = c("gradient boosted model",
                                             "linear model (all regressors)",
                                             "HINT best score seed 20",
                                             "Wellington best score seed 20"
                                             )
                                  )

# Plot the Matt CC Curve
png("../figures/seed20MCC.filtered.10M.png")
plot.mattcc.curve(all.stats.df) + theme_minimal(base_size = 15) + colScale
dev.off()

# Plot the ROC curve
png("../figures/seed20ROC.filtered.10M.png")
plot.roc.curve(all.stats.df, TRUE) + theme_minimal(base_size = 15) + colScale
dev.off()

# Plot the Prec-Rec curve
png("../figures/seed20PreRec.filtered.10M.png")
plot.precrecall.curve(all.stats.df) + theme_minimal(base_size = 15) + colScale
dev.off()
