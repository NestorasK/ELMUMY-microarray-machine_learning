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
)[3]
fsmeta_test <- list.files(
    path = path2read,
    pattern = ".*holdout.*progressing_MGUS.*",
    full.names = TRUE, recursive = TRUE
)[3]
path2save <- paste0(
    "results/experiments_caret/multiple_myeloma_progression/",
    "interpretation_optimizing_auc/variable_importance/"
)
path2read_models <- paste0(
    "results/experiments_caret/multiple_myeloma_progression/",
    "training_optimizing_auc/"
)

# Calculations ####
importance_all <- vector(mode = "list", length = length(fsmeta_train))
j <- 1
for (j in seq_len(length.out = length(fsmeta_train))) {
    cat("Use training metadata:", j, "out of", length(fsmeta_train), "\n")
    cat("Use training metadata:", fsmeta_train[j], "\n")

    fmeta_train <- fsmeta_train[j]
    meta_train <- fread(fmeta_train)

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

        # Variable importance
        import_transformi <- lapply(
            X = models_all,
            FUN = function(mli) {
                cat(
                    "Variable importance of", mli$method, "&",
                    transformationi, "...\n"
                )
                importi_list <- varImp(object = mli)
                pdf(
                    file = paste0(
                        path2save, mli$method, "_", transformationi,
                        ".pdf"
                    ), width = 5, height = 5
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

# Add ratio probes information
ratio_probes <- fread(
    paste0(
        "data/processed_gpl96_gpl570_affy44_platform/",
        "background_corrected_expression_ratios_onlyprobes.csv"
    )
)
colnames(ratio_probes) <- "ratio_probes"
ratio_probes$rn <- paste0("V", seq_len(length.out = nrow(ratio_probes)))
test <- merge.data.table(
    x = importance_all[, "rn"], y = ratio_probes,
    by = "rn", all = TRUE
)
test$rn_all <- ifelse(is.na(test$ratio_probes),
    yes = test$rn, no = test$ratio_probes
)
importance_all_plus_ratioprobes <- merge.data.table(
    x = unique(test[, c("rn", "rn_all")]),
    y = importance_all, by = "rn"
)
importance_all_plus_ratioprobes[, rn := NULL]
colnames(importance_all_plus_ratioprobes)[1] <- "rn"
fwrite(
    x = importance_all_plus_ratioprobes,
    file = paste0(path2save, "variable_importance_caret.csv")
)
