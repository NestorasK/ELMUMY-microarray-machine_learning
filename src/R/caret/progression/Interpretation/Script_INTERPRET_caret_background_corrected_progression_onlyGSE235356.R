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
path2save <- paste0(
    "results/experiments_caret/multiple_myeloma_progression_onlyGSE235356/",
    "interpretation_optimizing_auc/variable_importance/"
)
path2read_models <- paste0(
    "results/experiments_caret/multiple_myeloma_progression_onlyGSE235356/",
    "training_optimizing_auc/"
)

models <- c("glmnet", "rf", "gbm", "svmRadial", "svmLinear2")
transformations <- c("rma", "ranking", "qnorm", "binary_0.5", "ratios")

# Calculations ####
models_transf <- expand.grid(models = models, transformations = transformations)

load_models_helper <- function(linej) {
    modeli <- linej[1]
    transformationi <- linej[2]
    file_modeli <- paste0(
        path2read_models, modeli,
        "_", transformationi,
        "_FINAL.RData"
    )
    if (file.exists(file_modeli)) {
        cat("Loading:", file_modeli, "...\n")
        load(file_modeli)
        return(modeli_fit)
    } else {
        cat("No model found!\n")
    }
}
models_all <- apply(X = models_transf, MARGIN = 1, FUN = load_models_helper)
names(models_all) <- paste0(
    models_transf$models, "+",
    models_transf$transformations
)

# Variable importance
import_transformi <- lapply(
    X = seq_len(length.out = length(models_all)),
    FUN = function(j) {
        cat(
            "Variable importance of", names(models_all)[j], "...\n"
        )

        importi_list <- varImp(object = models_all[[j]])
        pdf(
            file = paste0(path2save, names(models_all)[j], ".pdf"),
            width = 5, height = 5
        )
        print(
            plot(importi_list,
                top = 20,
                main = names(models_all)[j]
            )
        )
        dev.off()
        importi <- data.table(importi_list$importance,
            keep.rownames = TRUE
        )
        importi[, c("method", "transformation") := tstrsplit(
            x = names(models_all)[j], split = "+", keep = c(1, 2),
            fixed = TRUE
        )]
        return(importi)
    }
)
importance_all <- rbindlist(l = import_transformi, fill = TRUE)

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
