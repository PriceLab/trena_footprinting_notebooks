# Load all libraries and functions

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

(.packages())

source("../my_R_functions/utility_functions.R")
source("../my_R_functions/stat_functions.R")
source("../my_R_functions/plot_functions.R")
source("/ssd/mrichard/github/BDDS/trenadb/src/utils.R")
source("/ssd/mrichard/github/BDDS/footprints/testdb/src/dbFunctions.R")

# Load the motif data and the big dataset, plus make the names better
load("../Rdata_files/motif_class_pairs.Rdata")
#load("/ssd/mrichard/data/joined.annotated.9.Rdata") # Takes ~2 min to load
all.TF.df.fimo.hint.well.annotated <- read.fst("/ssd/mrichard/data/joined.annotated.10M.Rdata")
colnames(all.TF.df.fimo.hint.well.annotated) <- make.names(colnames(all.TF.df.fimo.hint.well.annotated), unique=TRUE)

# Filter the entries w/o fp or hits
all.TF.df.fimo.hint.well.annotated %>%
    filter(h_max_score_16 > 3.6 | h_max_score_20 > 3.6 | w_min_score_16 < -2.3 | w_min_score_20 < -2.3) ->
    df_only_footprint_hits


# Remove columns we don't need for this part, then drop the original data
# Include the HINT/Wellington columns
cols_to_drop <- c('motifname', 'chrom', 'strand', 'loc',
                  'h_frac_16', 'h_frac_20', 'h_max_score_16', 'h_max_score_20',
                  'w_frac_16', 'w_frac_20', 'w_min_score_16', 'w_min_score_20')

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

remove(df_only_footprint_hits, all.TF.df.fimo.hint.well.annotated)

# Divide into train/test/validation sets and split X/Y
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

# MODEL 1: GBM
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
xgb.save(gbdt_medium, "../saved_models/xgboost_TF_site_predict_motif_only.10M.filtered.model")

# Create importance matrix plot
motif.class$class <- lapply(motif.class$class, make.names, unique=TRUE)
importance_matrix <- xgb.importance(colnames(X_train),model=gbdt_medium)
df <- as_data_frame(importance_matrix)
df.tf <- subset(df, Feature %in% unique(motif.class$class))
df.notf <- subset(df, !(Feature %in% unique(motif.class$class)))
tfclass.row <- c("TF_class", unname(as.list(colSums(df.tf[!(colnames(df.tf) %in% c("Feature"))]))) )
names(tfclass.row) <- colnames(df)
df.sum <- rbind(df.notf,tfclass.row)


png("../figures/motifOnlyImpMatrix.filtered.png")
ggplot(data=df.sum, aes(x=reorder(Feature, Gain), y=Gain)) +
    geom_bar(stat="identity") +
    coord_flip() +
    theme_minimal(base_size = 30) +
    labs(x = "Feature", y="Gain")
dev.off()

# Create stats for GBM
medium_pred_df <- make.pred.df.from.model(gbdt_medium, X_test, y_test)
colnames(medium_pred_df)[1] <- "ChIPseq.bound"
medium_stat_df <- make.stats.df.from.preds(medium_pred_df)


# MODEL 2: Full Linear Model
X_train_lin <- X_train
y_train_lin <- y_train
X_test_lin  <- X_test
y_test_lin  <- y_test

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

# Create stats for Full Linear Model
glm.pred.df <- make.pred.df.from.glm(glm.all, glm.df.test)
glm.stat.df <- make.stats.df.from.preds(glm.pred.df)


# MODEL 3: Individual Linear Models
stats.regressors.df <- data.frame()

for (this.regressor in colnames(X_train_lin)) {

    if (this.regressor %in% unique(motif.class$class)) {
        glm.formula <- paste("ChIPseq.bound ~ ", "as.factor(", this.regressor, ")",sep='')
    } else {
        glm.formula <- paste("ChIPseq.bound ~ ", this.regressor,sep='')
    }

    glm.df.train <- as.data.frame(cbind(y_train_lin, X_train_lin)) %>% rename("cs_hit" = "ChIPseq.bound")
    glm.df.test <-  as.data.frame(cbind(y_test_lin, X_test_lin)) %>% rename("cs_hit" = "ChIPseq.bound")

    # Make the model
    glm.single <- glm(as.formula(glm.formula), data=glm.df.train, family=binomial)
    glm.single$Model.Name <- paste("glm ", this.regressor, sep='')

    # Create the stats
    glm.pred.single.df <- make.pred.df.from.glm(glm.single, glm.df.test)
    glm.stat.single.df <- make.stats.df.from.preds(glm.pred.single.df)

    stats.regressors.df <- rbind(stats.regressors.df, glm.stat.single.df)

}

stats.regressors.df %>%
    filter(Model.Name %in% c("glm motifscore",
                             "glm gc_content",
                             "glm asinh_tss_dist")) ->
    stats.regressors.filtered.df

# Compare performance using 3 figures:
all.stats.df <- rbind(
    medium_stat_df,
    glm.stat.df,
    stats.regressors.filtered.df
    )

# Save the "all.stats.df" for later
save(all.stats.df, file = "/ssd/mrichard/data/motifOnlyAllStats.10M.filtered.Rdata")

# MCC curves
png("../figures/motifOnlyMCC.filtered.png")
plot.mattcc.curve(all.stats.df) + theme_minimal(base_size = 15)
dev.off()


# ROC curves
png("../figures/motifOnlyROC.filtered.png")
plot.roc.curve(all.stats.df) + theme_minimal(base_size = 15)
dev.off()


# Precision-Recall curves
png("../figures/motifOnlyPreRec.filtered.png")
plot.precrecall.curve(all.stats.df) + theme_minimal(base_size = 15)
dev.off()
