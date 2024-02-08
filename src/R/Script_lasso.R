rm(list = ls())
library(data.table)
library(glmnet)
library(pROC)

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
fmeta_train <- "data/processed_affy44_platform/metadata_train.csv"
fmeta_test <- "data/processed_affy44_platform/metadata_holdout.csv"
path2save <- "results/lasso_affy44_platform/"
reps <- 20

# Calculations ####
aucs_holdout <- matrix(data = NA, nrow = length(files_expr), ncol = reps)
aucs_mean_cv <- aucs_holdout
numfeatures <- aucs_holdout

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
            family = "binomial",
            type.measure = "auc"
        )
        aucs_mean_cv[ifile, repi] <- lasso_model$cvm[lasso_model$index[2]]
        numfeatures[ifile, repi] <- lasso_model$nzero[lasso_model$index[2]]

        # Test on hold out
        holdout_y_pred <- predict(
            object = lasso_model,
            newx = datain$holdout_x, type = "response"
        )
        aucs_holdout[ifile, repi] <- auc(
            response = datain$holdout_y, predictor = as.numeric(holdout_y_pred)
        )
    }
}
transformations <- sub(
    ".csv", "",
    sub(
        pattern = "expression_", replacement = "", x = basename(files_expr)
    )
)
aucs_mean_cv_dt <- data.table(aucs_mean_cv)
colnames(aucs_mean_cv_dt) <- paste0("rep_", 1:reps)
aucs_mean_cv_dt$transformations <- transformations
aucs_mean_cv_dt$metric <- "auc_cvmean"

aucs_holdout <- data.table(aucs_holdout)
colnames(aucs_holdout) <- paste0("rep_", 1:reps)
aucs_holdout$transformations <- transformations
aucs_holdout$metric <- "auc_holdout"

aucs_all <- rbindlist(l = list(aucs_mean_cv_dt, aucs_holdout))
fwrite(x = aucs_all, file = paste0(path2save, "aucs_all.csv"))

numfeatures <- data.table(numfeatures)
colnames(numfeatures) <- paste0("rep_", 1:reps)
numfeatures$transformations <- transformations
fwrite(x = numfeatures, file = paste0(path2save, "numfeatures.csv"))

# Plots
plot_perf <- ggplot(
    data = melt.data.table(
        data = aucs_all, id.vars = c("transformations", "metric"),
        value.name = "auc"
    ),
    mapping = aes(x = transformations, y = auc, fill = metric)
) +
    geom_boxplot() +
    ylim(0.5, 1) +
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
    plot = plot_numfeat, width = 7, height = 5
)
