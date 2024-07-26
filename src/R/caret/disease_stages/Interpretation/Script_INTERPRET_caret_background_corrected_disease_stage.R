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
# library(httpgd)
# hgd()
# hgd_browse()


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
fsmeta_train <- list.files(
    path = path2read, pattern = "metadata_train_classes", full.names = TRUE
)[5]
fsmeta_test <- list.files(
    path = path2read,
    pattern = ".*holdout.*progressing_MGUS.*",
    full.names = TRUE, recursive = TRUE
)[5]
path2save <- paste0(
    "results/experiments_caret/multiple_myeloma_stage/",
    "interpretation_optimizing_multiclass_auc/"
)
path2read_models <- paste0(
    "results/experiments_caret/multiple_myeloma_stage/",
    "training_optimizing_multiclass_auc/"
)

# Calculations ####
importance_all <- vector(mode = "list", length = length(fsmeta_train))
j <- 1
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
    import_all_transform <- foreach(
        ifile = seq_len(length.out = length(files_expr)), .combine = rbind
    ) %do% {
        fi <- files_expr[ifile]
        transformationi <- sub(
            ".csv", "",
            sub(
                pattern = "background_corrected_expression_",
                replacement = "", x = basename(fi)
            )
        )

        # Test models
        models <- c("glmnet", "rf", "gbm", "svmRadial", "svmLinear2")
        load_models_helper <- function(modeli) {
            file_modeli <- paste0(
                path2read_models, modeli, "_", transformationi,
                ".RData"
            )
            if (file.exists(file_modeli)) {
                cat("Loading", file_modeli, "...\n")
                load(file_modeli)

                return(modeli_fit)
            } else {
                cat("No model found!\n")
            }
        }
        models_all <- lapply(X = models, FUN = load_models_helper)
        names(models_all) <- models

        # Predictions all
        import_transformi <- lapply(
            X = models_all,
            FUN = function(mli) {
                cat(
                    "Variable importance of", mli$method, "&",
                    transformationi, "...\n"
                )
                importi_list <- varImp(object = mli)
                if (mli$method %in% c("gbm", "rf")) {
                    widthi <- 5
                } else {
                    widthi <- 10
                }
                pdf(
                    file = paste0(
                        path2save, mli$method, "_", transformationi,
                        ".pdf"
                    ), width = widthi, height = 5
                )
                print(
                    plot(importi_list,
                        top = 20,
                        main = paste(mli$method, transformationi, sep = " | ")
                    )
                )
                dev.off()
                importi <- data.table(importi_list$importance,
                    keep.rownames = TRUE
                )
                importi$method <- mli$method
                return(importi)
            }
        )
        import_transformi_dt <- rbindlist(l = import_transformi, fill = TRUE)
        import_transformi_dt$transformation <- transformationi
        return(import_transformi_dt)
    }
    import_all_transform$meta_train <- basename(fmeta_train)
    import_all_transform$transformation[
        grepl(pattern = "qnorm_", x = import_all_transform$transformation)
    ] <- "qnorm"
    importance_all[[j]] <- import_all_transform
}
importance_all <- rbindlist(l = importance_all)
fwrite(
    x = importance_all,
    file = paste0(path2save, "variable_importance_caret.csv")
)
