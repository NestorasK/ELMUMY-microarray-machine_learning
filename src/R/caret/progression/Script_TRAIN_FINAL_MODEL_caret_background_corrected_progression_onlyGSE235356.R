rm(list = ls())
library(data.table)
library(dplyr)
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
files_expr_stable <- paste0(
    path2read,
    c(
        "background_corrected_expression_rma.csv",
        # "background_corrected_expression_binary_0.csv",
        # "background_corrected_expression_binary_0.25.csv",
        "background_corrected_expression_binary_0.5.csv",
        # "background_corrected_expression_binary_0.75.csv",
        "background_corrected_expression_ranking.csv",
        "background_corrected_expression_ratios.csv"
    )
)
metadata_all <- fread(
    paste0(
        "data/processed_gpl96_gpl570_affy44_platform/",
        "metadata_holdout_classes:['MGUS', 'MM', 'progressing_MGUS']_",
        "dataset:['EMTAB317'].csv"
    )
)
path2save <- paste0(
    "results/experiments_caret/",
    "multiple_myeloma_progression_onlyGSE235356/"
)

# Models
models <- c("glmnet", "rf", "gbm", "svmLinear2", "svmRadial")

# Parameters tuning
fit_control <- trainControl( ## 10-fold CV
    method = "repeatedcv",
    number = 10,
    ## repeated ten times
    repeats = 10,
    classProbs = TRUE,
    summaryFunction = twoClassSummary
)

# Calculations ####
metadata_gse235356 <- metadata_all[dataset == "GSE235356", ]
registerDoMC(cores = detectCores())
set.seed(42)
for (j in seq_len(length.out = length(files_expr_stable))) {
    cat("\n###########################################")
    cat("\nTraining using:", files_expr_stable[j], "\n")
    cat("###########################################\n")
    expr <- fread(files_expr_stable[j])
    transformationj <- sub(
        ".csv", "",
        sub(
            pattern = "background_corrected_expression_",
            replacement = "", x = basename(files_expr_stable[j])
        )
    )
    if (transformationj == "ratios") {
        expr_t <- transpose(l = expr, keep.names = "rn")
    } else {
        expr_t <- transpose(l = expr, keep.names = "rn", make.names = "rn")
    }
    expr_t_meta <- merge.data.table(
        x = metadata_gse235356[, -"dataset"],
        y = expr_t, by = "rn"
    )
    train_x <- as.matrix(x = expr_t_meta[, -c("rn", "class")])
    train_y <- expr_t_meta[, class]

    # Train models
    train_helper <- function(modeli) {
        cat("\nTraining:", modeli, "\n")
        file_modeli <- paste0(
            path2save,
            "training_optimizing_auc/",
            modeli,
            "_", transformationj,
            "_FINAL.RData"
        )
        if (file.exists(file_modeli)) {
            cat("Loading model...\n")
            load(file_modeli)
            return(modeli_fit)
        }
        if (modeli == "svmRadial") {
            modeli_fit <- caret::train(
                x = train_x, y = train_y,
                method = "svmRadial",
                trControl = fit_control,
                tuneGrid = expand.grid(
                    sigma = kernlab::sigest(
                        train_x,
                        na.action = na.omit, scaled = TRUE
                    ),
                    C = c(0.25, 0.5, 2^c(1:10))
                ),
                metric = "ROC",
                verbose = FALSE
            )
        } else {
            modeli_fit <- caret::train(
                x = train_x, y = train_y,
                method = modeli,
                preProcess = c("center", "scale"),
                trControl = fit_control,
                tuneLength = 10,
                metric = "ROC",
                verbose = FALSE
            )
        }
        cat("Saving:", file_modeli, "\n")
        save(modeli_fit, file = file_modeli)
        pdf(
            file = paste0(
                path2save, "training_optimizing_auc/",
                modeli, "_tuning_",
                transformationj,
                "_FINAL.pdf"
            ),
            width = 10, height = 5
        )
        print(
            plot(
                modeli_fit,
                main = paste(modeli, transformationj, sep = " | ")
            )
        )
        dev.off()
        return(modeli_fit)
    }
    models_all <- lapply(X = models, FUN = train_helper)
    names(models_all) <- models

    # Compare all
    resamps <- resamples(models_all)
    models_comparison <- summary(resamps)
    save(models_comparison,
        file = paste0(
            path2save,
            "training_optimizing_auc/",
            "models_comparison_",
            transformationj,
            "_FINAL.RData"
        )
    )
    pdf(
        file = paste0(
            path2save, "training_optimizing_auc/",
            "models_comparison_",
            transformationj,
            "_FINAL.pdf"
        ),
        width = 8,
        height = 4
    )
    theme1 <- trellis.par.get()
    theme1$plot.symbol$col <- rgb(.2, .2, .2, .4)
    theme1$plot.symbol$pch <- 16
    theme1$plot.line$col <- rgb(1, 0, 0, .7)
    theme1$plot.line$lwd <- 2
    trellis.par.set(theme1)
    print(
        bwplot(
            resamps,
            layout = c(3, 1), main = "Cross Validation performance"
        )
    )
    dev.off()
}
print(warnings())
