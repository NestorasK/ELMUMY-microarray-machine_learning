rm(list = ls())
library(data.table)
library(preprocessCore)
library(doMC)
library(glmnet)
library(caret)
library(pROC)
library(ggplot2)
library(httpgd)
# hgd()
# hgd_browse()

# Input ####
path2read <- "data/processed_gpl96_gpl570_affy44_platform/"
file_expr_stable <- paste0(
    path2read,
    "background_corrected_expression_rma.csv"
)
metadata_all <- fread(
    "data/processed_gpl96_gpl570_affy44_platform/metadata_holdout_classes:['MGUS', 'MM', 'progressing_MGUS']_dataset:['EMTAB317'].csv"
)
path2save <- "results/processed_gpl96_gpl570_affy44_platform/"
reps <- 20

# Calculations ####
metadata_gse235356 <- metadata_all[dataset == "GSE235356", ]
set.seed(42)

cat("\nTraining using:", file_expr_stable, "\n")
expr <- fread(file_expr_stable)
expr_t <- transpose(l = expr, keep.names = "rn", make.names = "rn")
expr_t_meta <- merge.data.table(
    x = metadata_gse235356[, -"dataset"],
    y = expr_t, by = "rn"
)

registerDoMC(cores = reps)
model_out <- foreach(repi = seq_len(reps), .combine = rbind) %dopar% {
    # make holdout and train folds
    ncvfolds <- caret::createDataPartition(
        y = expr_t_meta[, class],
        p = 0.7
    )
    holdout_x <- as.matrix(
        x = expr_t_meta[-unlist(ncvfolds), -c("rn", "class")]
    )
    rownames(holdout_x) <- expr_t_meta[-unlist(ncvfolds), rn]
    holdout_y <- expr_t_meta[-unlist(ncvfolds), class]

    train_x <- as.matrix(
        x = expr_t_meta[unlist(ncvfolds), -c("rn", "class")]
    )
    rownames(train_x) <- expr_t_meta[unlist(ncvfolds), rn]
    train_y <- expr_t_meta[unlist(ncvfolds), class]

    # Quantile normalize
    train_x_qnorm <- t(normalize.quantiles(t(train_x)))
    rownames(train_x_qnorm) <- rownames(train_x)
    colnames(train_x_qnorm) <- colnames(train_x)
    holdout_x_qnorm <- t(
        x = normalize.quantiles.use.target(
            x = t(holdout_x),
            target = train_x_qnorm[1, ]
        )
    )
    colnames(holdout_x_qnorm) <- colnames(holdout_x)
    rownames(holdout_x_qnorm) <- rownames(holdout_x)

    # Stratified CV folds
    folds_repi <- caret::createFolds(y = train_y, list = FALSE)

    # train lasso
    lasso_model <- cv.glmnet(
        x = train_x, y = train_y,
        type.measure = "auc", foldid = folds_repi, family = "binomial"
    )
    auc_cv <- lasso_model$cvm[lasso_model$index[2]]
    number_of_features <- lasso_model$nzero[lasso_model$index[2]]
    names(number_of_features) <- NULL

    # Test on hold out
    holdout_y_prob <- predict(
        object = lasso_model,
        newx = holdout_x, type = "response"
    )
    auc_holdout <- pROC::auc(
        response = holdout_y,
        predictor = holdout_y_prob[, 1]
    )
    return(
        c(
            "number_of_features" = number_of_features,
            "auc_cv" = auc_cv, "auc_holdout" = auc_holdout
        )
    )
}
model_out_dt <- data.table(model_out, keep.rownames = TRUE)
model_out_dt$transformations <- "qnorm"
filename <- paste0(path2save, "GSE235356_alone_lasso_output_qnorm.csv")
fwrite(x = model_out_dt, file = filename)

# Read all other models but qnorm
models_output <- fread(
    paste0(path2save, "GSE235356_alone_lasso_output.csv")
)
models_outputsdt <- rbindlist(l = list(models_output, model_out_dt))

# Plots
plot_perf <- ggplot(
    data = melt.data.table(
        data = models_outputsdt[
            , c("auc_cv", "auc_holdout", "transformations")
        ],
        id.vars = "transformations", variable.name = "metric",
        value.name = "auc"
    ),
    mapping = aes(x = transformations, y = auc, fill = metric)
) +
    geom_boxplot() +
    ggtitle("GSE235356 - MGUS vs progressive MGUS")

ggsave(
    filename = paste0(path2save, "GSE235356_alone_lasso_perf.pdf"),
    plot = plot_perf, height = 5, width = 10
)

plot_numfeat <- ggplot(
    data = models_outputsdt[, c("transformations", "number_of_features")],
    mapping = aes(x = transformations, y = number_of_features)
) +
    geom_boxplot() +
    ggtitle("GSE235356 - MGUS vs progressive MGUS")

ggsave(
    filename = paste0(path2save, "GSE235356_alone_lasso_numfeatures.pdf"),
    plot = plot_numfeat, height = 5, width = 6
)
