rm(list = ls())
library(data.table)
library(glmnet)
library(pROC)
library(ggplot2)
library(httpgd)
hgd()
hgd_browse()

# Input ####
path2read <- "data/processed_gpl96_gpl570_affy44_platform"
files_expr <- list.files(
    path = path2read, pattern = "background_corrected",
    full.names = TRUE
)
fmeta_train <- "data/processed_gpl96_gpl570_affy44_platform/metadata_train.csv"
fmeta_test <- "data/processed_gpl96_gpl570_affy44_platform/metadata_holdout.csv"
path2save <- "results/processed_gpl96_gpl570_affy44_platform/"
reps <- 20

# Calculations ####
meta_train <- fread(fmeta_train)
meta_holdout <- fread(fmeta_test)

set.seed(42)
accuracys_mean_cv <- matrix(data = NA, nrow = length(files_expr), ncol = reps)
numfeatures <- accuracys_mean_cv
accuracys_holdout <- vector("list", length = length(files_expr))
for (ifile in seq_len(length.out = length(files_expr))) {
    fi <- files_expr[ifile]
    cat("\nTraining using:", fi, "\n")
    expr <- fread(fi)

    source("src/R/fetch_data.R")
    datain <- fetch_data(
        expres = expr,
        meta_train = meta_train, meta_holdout = meta_holdout
    )
    transformationi <- sub(
        ".csv", "",
        sub(
            pattern = "background_corrected_expression_",
            replacement = "", x = basename(fi)
        )
    )
    accuracys_holdout_repi <- vector("list", reps)
    for (repi in 1:reps) {
        cat("-", repi, "out of", reps, "\n")

        # fit lasso
        lasso_model <- cv.glmnet(
            x = datain$train_x, y = datain$train_y,
            family = "multinomial",
            type.measure = "class"
        )
        # Plot cv error
        # plot(lasso_model, main = paste0(transformationi, "\n"))
        accuracys_mean_cv[ifile, repi] <- (
            1 - lasso_model$cvm[lasso_model$index[2]]
        )
        numfeatures[ifile, repi] <- lasso_model$nzero[lasso_model$index[2]]

        # Test on hold out
        holdout_y_pred <- predict(
            object = lasso_model,
            newx = datain$holdout_x, type = "class"
        )
        rownames(holdout_y_pred) <- rownames(datain$holdout_x)
        preds_dt <- merge.data.table(
            x = meta_holdout,
            y = data.table(holdout_y_pred, keep.rownames = TRUE),
            by = "rn"
        )
        colnames(preds_dt)[4] <- "class_hat"
        accuracy_holdout_per_dataset <- preds_dt[
            , .(accuracy = mean(class_hat == class)),
            by = dataset
        ]
        accuracys_holdout_all <- rbindlist(
            list(
                data.table(
                    dataset = "all",
                    accuracy = mean(preds_dt$class_hat == preds_dt$class)
                ),
                accuracy_holdout_per_dataset
            )
        )
        accuracys_holdout_all$transformations <- transformationi
        accuracys_holdout_repi[[repi]] <- accuracys_holdout_all
    }
    accuracys_holdout[[ifile]] <- rbindlist(l = accuracys_holdout_repi)
}
accuracys_holdout <- rbindlist(l = accuracys_holdout)
accuracys_holdout$metric <- "accuracy_holdout"

# accuracy mean CV
transformations <- sub(
    ".csv", "",
    sub(
        pattern = "background_corrected_expression_",
        replacement = "", x = basename(files_expr)
    )
)
accuracys_mean_cv_dt <- data.table(accuracys_mean_cv)
colnames(accuracys_mean_cv_dt) <- paste0("rep_", 1:reps)
accuracys_mean_cv_dt$transformations <- transformations
accuracys_mean_cv_dt$metric <- "accuracy_cvmean"
accuracys_mean_cv_dt <- melt.data.table(
    data = accuracys_mean_cv_dt,
    id.vars = c("metric", "transformations"), value.name = "accuracy"
)
accuracys_mean_cv_dt$dataset <- "all"
accuracys_all <- rbindlist(
    l = list(
        accuracys_holdout, accuracys_mean_cv_dt[, -"variable"]
    ),
    use.names = TRUE
)
fwrite(
    x = accuracys_all,
    file = paste0(path2save, "background_corrected_accuracys_all.csv")
)

numfeatures <- data.table(numfeatures)
colnames(numfeatures) <- paste0("rep_", 1:reps)
numfeatures$transformations <- transformations
fwrite(
    x = numfeatures,
    file = paste0(path2save, "background_corrected_numfeatures.csv")
)

# Plots
plot_perf <- ggplot(
    data = accuracys_all,
    mapping = aes(x = transformations, y = accuracy, fill = metric)
) +
    geom_boxplot() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom"
    ) +
    facet_grid(cols = vars(dataset)) +
    scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1))
ggsave(
    filename = paste0(
        path2save, "background_corrected_boxplot_performance.pdf"
    ),
    plot = plot_perf, width = 15, height = 5
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
    filename = paste0(
        path2save, "background_corrected_boxplot_numfeatures.pdf"
    ),
    plot = plot_numfeat, width = 4, height = 4
)
