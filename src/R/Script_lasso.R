rm(list = ls())
library(data.table)
library(glmnet)
library(pROC)
library(ggplot2)

# Input ####
files_expr <- c(
    "data/processed_gpl96_gpl570_affy44_platform/expression_rma.csv",
    "data/processed_gpl96_gpl570_affy44_platform/expression_binary_0.csv",
    "data/processed_gpl96_gpl570_affy44_platform/expression_binary_0.25.csv",
    "data/processed_gpl96_gpl570_affy44_platform/expression_binary_0.5.csv",
    "data/processed_gpl96_gpl570_affy44_platform/expression_binary_0.75.csv",
    "data/processed_gpl96_gpl570_affy44_platform/expression_ranking.csv",
    "data/processed_gpl96_gpl570_affy44_platform/expression_ratios.csv",
    "data/processed_gpl96_gpl570_affy44_platform/expression_ratiosfromranks.csv"
)
fmeta_train <- "data/processed_glp96_gpl570_platform/metadata_train.csv"
fmeta_test <- "data/processed_glp96_gpl570_platform/metadata_holdout.csv"
path2save <- "results/processed_glp96_gpl570_platform/"
reps <- 20

# Calculations ####
accuracys_holdout <- matrix(data = NA, nrow = length(files_expr), ncol = reps)
accuracys_mean_cv <- accuracys_holdout
numfeatures <- accuracys_holdout

for (ifile in seq_len(length.out = length(files_expr))) {
    fi <- files_expr[ifile]
    cat("\nTraining using", fi, "\n")
    expr <- fread(fi)
    meta_train <- fread(fmeta_train)
    meta_holdout <- fread(fmeta_test)

    source("src/R/fetch_data.R")
    datain <- fetch_data(
        expres = expr,
        meta_train = meta_train, meta_holdout = meta_holdout
    )
    for (repi in 1:reps) {
        cat("-", repi, "out of", reps, "\n")

        # fit lasso
        lasso_model <- cv.glmnet(
            x = datain$train_x, y = datain$train_y,
            family = "multinomial",
            type.measure = "class"
        )
        accuracys_mean_cv[ifile, repi] <- (
            1 - lasso_model$cvm[lasso_model$index[2]]
        )
        numfeatures[ifile, repi] <- lasso_model$nzero[lasso_model$index[2]]

        # Test on hold out
        holdout_y_pred <- predict(
            object = lasso_model,
            newx = datain$holdout_x, type = "class"
        )
        accuracys_holdout[ifile, repi] <- mean(
            as.vector(holdout_y_pred) == datain$holdout_y
        )
    }
}
transformations <- sub(
    ".csv", "",
    sub(
        pattern = "expression_", replacement = "", x = basename(files_expr)
    )
)
accuracys_mean_cv_dt <- data.table(accuracys_mean_cv)
colnames(accuracys_mean_cv_dt) <- paste0("rep_", 1:reps)
accuracys_mean_cv_dt$transformations <- transformations
accuracys_mean_cv_dt$metric <- "accuracy_cvmean"

accuracys_holdout <- data.table(accuracys_holdout)
colnames(accuracys_holdout) <- paste0("rep_", 1:reps)
accuracys_holdout$transformations <- transformations
accuracys_holdout$metric <- "accuracy_holdout"

accuracys_all <- rbindlist(l = list(accuracys_mean_cv_dt, accuracys_holdout))
fwrite(x = accuracys_all, file = paste0(path2save, "accuracys_all.csv"))

accuracys_all <- fread(paste0(path2save, "accuracys_all.csv"))


numfeatures <- data.table(numfeatures)
colnames(numfeatures) <- paste0("rep_", 1:reps)
numfeatures$transformations <- transformations
fwrite(x = numfeatures, file = paste0(path2save, "numfeatures.csv"))

numfeatures <- fread(paste0(path2save, "numfeatures.csv"))

# Plots
plot_perf <- ggplot(
    data = melt.data.table(
        data = accuracys_all, id.vars = c("transformations", "metric"),
        value.name = "accuracy"
    ),
    mapping = aes(x = transformations, y = accuracy, fill = metric)
) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(
    filename = paste0(path2save, "boxplot_performance.pdf"),
    plot = plot_perf, width = 7, height = 5
)

plot_numfeat <- ggplot(
    data = melt.data.table(
        data = numfeatures, id.vars = c("transformations"),
        value.name = "num_features"
    ),
    mapping = aes(x = transformations, y = num_features)
) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(
    filename = paste0(path2save, "boxplot_numfeatures.pdf"),
    plot = plot_numfeat, width = 4, height = 4
)
