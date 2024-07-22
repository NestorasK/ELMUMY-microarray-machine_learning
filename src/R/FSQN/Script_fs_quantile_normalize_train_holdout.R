# Based on https://pubmed.ncbi.nlm.nih.gov/29360996/
# DO NOT USE !!!!
# DIFFERENT NUMBER OF SAMPLES GIVE DIFFERENT VALUES PER SAMPLE !!!!
rm(list = ls())
library(data.table)
library(preprocessCore)
library(ggplot2)
library(httpgd)
hgd()
hgd_browse()
source("src/R/FSQN/fsqn.R")

# Input ####
path2read <- "data/processed_gpl96_gpl570_affy44_platform/"
file_expr <- paste0(
    path2read,
    "background_corrected_expression_rma.csv"
)
fsmeta_train <- list.files(
    path = path2read, pattern = "metadata_train.*Normal*",
    full.names = TRUE, recursive = TRUE
)
cat("Metadata to train:", fsmeta_train, "\n")
fsmeta_test <- list.files(
    path = path2read,
    pattern = ".*holdout.*Normal.*",
    full.names = TRUE, recursive = TRUE
)
cat("Metadata to test:", fsmeta_test, "\n")

# Calculations ####
expr <- fread(file_expr)

meta_train <- fread(fsmeta_train[1])
meta_holdout <- fread(fsmeta_test[1])

source("src/R/fetch_data.R")
datain <- fetch_data(
    expres = expr,
    meta_train = meta_train, meta_holdout = meta_holdout
)

# train
train_x_qnorm <- datain$train_x
cat("\nTraining data rma\n")
print(train_x_qnorm[1:5, 1:5])

# holdout



# test it only on 5 genes
holdout_x_fsqnorm_2samples <- quantileNormalizeByFeature(
    matrix_to_normalize = datain$holdout_x[1:2, 1:5],
    target_distribution_matrix = datain$train_x[, 1:5]
)
holdout_x_fsqnorm_2samples


holdout_x_fsqnorm_5samples <- quantileNormalizeByFeature(
    matrix_to_normalize = datain$holdout_x[1:5, 1:5],
    target_distribution_matrix = datain$train_x[, 1:5]
)
holdout_x_fsqnorm_5samples


holdout_x_fsqnorm_10samples <- quantileNormalizeByFeature(
    matrix_to_normalize = datain$holdout_x[1:10, 1:5],
    target_distribution_matrix = datain$train_x[, 1:5]
)
holdout_x_fsqnorm_10samples
