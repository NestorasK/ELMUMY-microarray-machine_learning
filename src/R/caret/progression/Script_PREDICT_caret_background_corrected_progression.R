rm(list = ls())
library(data.table)
library(caret)
library(randomForest)
library(glmnet)
library(gbm)
library(xgboost)
library(e1071)
library(kernlab)
library(doMC)
library(pROC)

# Input ####
path2read <- "data/processed_gpl96_gpl570_affy44_platform/"
files_expr_stable <- paste0(
    path2read,
    c(
        "background_corrected_expression_rma.csv",
        "background_corrected_expression_binary_0.csv",
        # "background_corrected_expression_binary_0.25.csv",
        "background_corrected_expression_binary_0.5.csv",
        # "background_corrected_expression_binary_0.75.csv",
        "background_corrected_expression_ranking.csv",
        "background_corrected_expression_ratios.csv"
    )
)
fsmeta_train <- sort(
    list.files(
        path = path2read, pattern = "metadata_train_classes", full.names = TRUE
    )
)[-5]
fsmeta_test <- sort(
    list.files(
        path = path2read,
        pattern = ".*holdout.*progressing_MGUS.*",
        full.names = TRUE, recursive = TRUE
    )
)[-5]
path2save <- paste0(
    "results/experiments_caret/multiple_myeloma_progression/",
    "predictions_optimizing_auc/"
)
path2read_models <- paste0(
    "results/experiments_caret/multiple_myeloma_progression/",
    "training_optimizing_auc/"
)

# Calculations ####
preds_all <- vector(mode = "list", length = length(fsmeta_train))
for (j in seq_len(length.out = length(fsmeta_train))) {
    cat("\n\nUse training metadata:", j, "out of", length(fsmeta_train), "\n")
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
    )[1]
    files_expr <- c(files_expr_stable, files_expr_train)
    preds_all_transform <- foreach(
        ifile = seq_len(length.out = length(files_expr)), .combine = rbind
    ) %do% {
        fi <- files_expr[ifile]
        cat("\nTesting using:", fi, "\n\n")
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
        file_extension <- sub(
            pattern = ".csv", replacement = "",
            x = sub(
                pattern = "metadata_train_classes:",
                replacement = "", x = basename(fmeta_train)
            )
        )

        # Test models
        models <- c("glmnet", "rf", "gbm", "svmRadial", "svmLinear2")
        load_models_helper <- function(modeli) {
            file_modeli <- paste0(
                path2read_models, modeli,
                "_", transformationi,
                "_", file_extension,
                ".RData"
            )
            if (file.exists(file_modeli)) {
                cat("Loading:", file_modeli, "...\n")
                load(file_modeli)
                return(modeli_fit)
            } else {
                cat("No model found!\n")
            }
        }
        models_all <- lapply(X = models, FUN = load_models_helper)
        names(models_all) <- models

        # Predictions all
        preds_transformi <- lapply(
            X = models_all,
            FUN = function(mli) {
                cat("Predict", mli$method, "...\n")
                prd_mi <- data.frame(
                    predict.train(
                        object = mli, newdata = datain$holdout_x,
                        type = "prob"
                    )
                )
                prd_mi$method <- mli$method
                prd_mi$rn <- rownames(datain$holdout_x)
                return(prd_mi)
            }
        )
        preds_transformi_dt <- rbindlist(l = preds_transformi)
        preds_transformi_dt$transformation <- transformationi
        preds_transformi_dt$class_pred <- c("MGUS", "MM")[
            apply(
                X = as.matrix(
                    x = preds_transformi_dt[, c("MGUS", "MM")]
                ),
                MARGIN = 1,
                which.max
            )
        ]
        return(preds_transformi_dt)
    }
    preds_all_transform$meta_train <- basename(fmeta_train)
    preds_all_transform_out <- merge.data.table(
        x = meta_holdout, y = preds_all_transform, by = "rn"
    )
    preds_all_transform_out$transformation[
        grepl(pattern = "qnorm_", x = preds_all_transform_out$transformation)
    ] <- "qnorm"
    preds_all[[j]] <- preds_all_transform_out
}
preds_all <- rbindlist(l = preds_all)
fwrite(x = preds_all, file = paste0(path2save, "predictions.csv"))
