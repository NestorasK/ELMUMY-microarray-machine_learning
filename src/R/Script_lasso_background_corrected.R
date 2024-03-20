rm(list = ls())
library(data.table)
library(glmnet)
library(pROC)
library(ggplot2)
library(httpgd)
# hgd()
# hgd_browse()

# Input ####
path2read <- "data/processed_gpl96_gpl570_affy44_platform/"
files_expr_stable <- paste0(
    path2read,
    c(
        "background_corrected_expression_rma.csv",
        "background_corrected_expression_binary_0.csv",
        "background_corrected_expression_binary_0.25.csv",
        "background_corrected_expression_binary_0.5.csv",
        "background_corrected_expression_binary_0.75.csv",
        "background_corrected_expression_ranking.csv"
    )
)
fsmeta_train <- list.files(
    path = path2read, pattern = "metadata_train_classes", full.names = TRUE
)
fsmeta_test <- sub(
    pattern = "_train_", replacement = "_holdout_",
    x = fsmeta_train
)
path2save <- "results/processed_gpl96_gpl570_affy44_platform/"
reps <- 20

# Calculations ####
set.seed(42)
for (j in seq_len(length.out = length(fsmeta_train))) {
    cat("Use training metadata:", j, "out of", length(fsmeta_train), "\n")
    cat("Use training metadata:", fsmeta_train[j], "\n")

    fmeta_train <- fsmeta_train[j]
    fmeta_test <- fsmeta_test[j]

    meta_train <- fread(fmeta_train)
    meta_holdout <- fread(fmeta_test)

    files_expr_train <- grep(
        pattern = "_expression_",
        x = grep(
            pattern = unlist(
                strsplit(
                    x = fmeta_train, split = "classes:", fixed = TRUE
                )
            )[2],
            x = list.files(
                path = path2read,
                full.names = TRUE
            ),
            value = TRUE, fixed = TRUE
        ), value = TRUE, fixed = TRUE
    )
    files_expr <- c(files_expr_stable, files_expr_train)
    accuracys_mean_cv <- matrix(
        data = NA, nrow = length(files_expr),
        ncol = reps
    )
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
                family = "binomial",
                type.measure = "auc"
            )
            # Plot cv error
            # plot(lasso_model, main = paste0(transformationi, "\n"))
            accuracys_mean_cv[ifile, repi] <- (
                lasso_model$cvm[lasso_model$index[2]]
            )
            numfeatures[ifile, repi] <- lasso_model$nzero[lasso_model$index[2]]

            # Test on hold out
            holdout_y_class <- predict(
                object = lasso_model,
                newx = datain$holdout_x, type = "class"
            )
            holdout_y_prob <- predict(
                object = lasso_model,
                newx = datain$holdout_x, type = "response"
            )
            holdout_y_pred <- merge.data.table(
                x = data.table(holdout_y_class, keep.rownames = TRUE),
                y = data.table(holdout_y_prob, keep.rownames = TRUE),
                by = "rn"
            )
            colnames(holdout_y_pred)[c(2, 3)] <- c(
                "predicted_class", "predicted_prob"
            )
            preds_dt <- merge.data.table(
                x = meta_holdout,
                y = holdout_y_pred,
                by = "rn"
            )

            # Accuracy
            accuracy_holdout_per_dataset <- preds_dt[
                , .(performance = mean(class == predicted_class)),
                by = dataset
            ]
            accuracys_holdout_all <- rbindlist(
                list(
                    data.table(
                        dataset = "all",
                        performance = mean(
                            preds_dt$predicted_class == preds_dt$class
                        )
                    ),
                    accuracy_holdout_per_dataset
                )
            )
            accuracys_holdout_all$metric <- "accuracy"

            # AUCs
            auc_holdout_per_dataset <- preds_dt[
                !dataset %in% c("GSE235356", "GSE5900"),
                .(performance = pROC::auc(
                    response = class,
                    predictor = predicted_prob
                )),
                by = dataset
            ]
            aucs_holdout_all <- rbindlist(
                list(
                    data.table(
                        dataset = "all",
                        performance = pROC::auc(
                            response = preds_dt[
                                !dataset %in% c("GSE235356", "GSE5900"), class
                            ],
                            predictor = preds_dt[
                                !dataset %in% c("GSE235356", "GSE5900"),
                                predicted_prob
                            ]
                        )
                    ),
                    auc_holdout_per_dataset
                )
            )
            aucs_holdout_all$metric <- "AUC"

            # Combine all performance metrics
            aucs_holdout_all$performance <- as.numeric(
                aucs_holdout_all$performance
            )
            accuracys_holdout_all <- rbindlist(
                l = list(accuracys_holdout_all, aucs_holdout_all)
            )
            accuracys_holdout_all$transformations <- transformationi
            accuracys_holdout_repi[[repi]] <- accuracys_holdout_all
        }
        accuracys_holdout[[ifile]] <- rbindlist(l = accuracys_holdout_repi)
    }
    accuracys_holdout <- rbindlist(l = accuracys_holdout)
    accuracys_holdout$transformations[
        grepl(pattern = "qnorm", x = accuracys_holdout$transformations)
    ] <- "qnorm"
    accuracys_holdout$transformations[
        grepl(pattern = "ratios", x = accuracys_holdout$transformations)
    ] <- "ratios"

    # accuracy mean CV
    transformations <- unique(accuracys_holdout$transformations)
    accuracys_mean_cv_dt <- data.table(accuracys_mean_cv)
    colnames(accuracys_mean_cv_dt) <- paste0("rep_", 1:reps)
    accuracys_mean_cv_dt$transformations <- transformations
    accuracys_mean_cv_dt$metric <- "auc_cvmean"
    accuracys_mean_cv_dt <- melt.data.table(
        data = accuracys_mean_cv_dt,
        id.vars = c("metric", "transformations"), value.name = "accuracy"
    )
    accuracys_mean_cv_dt$dataset <- unlist(
        strsplit(
            x = unlist(
                strsplit(
                    x = basename(fmeta_train), split = "dataset:['",
                    fixed = TRUE
                )
            )[2],
            split = "'].csv", fixed = TRUE
        )
    )
    colnames(accuracys_mean_cv_dt)[4] <- "performance"

    # Combine performances
    accuracys_all <- rbindlist(
        l = list(
            accuracys_holdout, accuracys_mean_cv_dt[, -"variable"]
        ),
        use.names = TRUE
    )
    filename_perf <- paste0(
        path2save,
        sub(
            pattern = "metadata_train_classes",
            replacement = "performance_all", x = basename(fmeta_train)
        )
    )
    cat("Saving performance file: '", filename_perf, "'\n", sep = "")
    fwrite(x = accuracys_all, file = filename_perf)

    # number of features
    numfeatures <- data.table(numfeatures)
    colnames(numfeatures) <- paste0("rep_", 1:reps)
    numfeatures$transformations <- transformations
    filename_numfeatures <- sub(
        pattern = "performance_all",
        replacement = "numfeatures", x = filename_perf
    )
    cat("Saving numfeatures file: '", filename_numfeatures, "'\n", sep = "")
    fwrite(x = numfeatures, file = filename_numfeatures)

    # Plots ####
    # - performance
    plot_perf <- ggplot(
        data = accuracys_all[metric %in% c("AUC", "auc_cvmean"), ],
        mapping = aes(x = transformations, y = performance, fill = metric)
    ) +
        geom_boxplot() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "bottom"
        ) +
        facet_grid(cols = vars(dataset)) +
        scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1))

    filename_plotperf <- sub(
        pattern = "performance_all",
        replacement = "boxplot_performance_AUC", x = filename_perf
    )
    filename_plotperf <- sub(
        pattern = "csv",
        replacement = "pdf", x = filename_plotperf
    )
    ggsave(
        filename = filename_plotperf, plot = plot_perf, width = 15,
        height = 5
    )

    # - number of features
    plot_numfeat <- ggplot(
        data = melt.data.table(
            data = numfeatures, id.vars = c("transformations"),
            value.name = "num_features"
        ),
        mapping = aes(x = transformations, y = num_features)
    ) +
        geom_boxplot() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

    filename_plotnumfeat <- sub(
        pattern = "boxplot_performance",
        replacement = "boxplot_numfeatures", x = filename_plotperf
    )
    ggsave(
        filename = filename_plotnumfeat,
        plot = plot_numfeat, width = 4, height = 4
    )
}
