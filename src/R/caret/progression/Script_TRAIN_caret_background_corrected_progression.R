rm(list = ls())
library(data.table)
library(caret)
library(randomForest)
library(glmnet)
library(gbm)
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
fsmeta_train <- list.files(
    path = path2read, pattern = "metadata_train_classes", full.names = TRUE
)[3]
fsmeta_test <- list.files(
    path = path2read,
    pattern = ".*holdout.*progressing_MGUS.*",
    full.names = TRUE, recursive = TRUE
)[4]
path2save <- paste0(
    "results/experiments_caret/",
    "multiple_myeloma_progression/training_optimizing_auc/"
)
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
doMC::registerDoMC(cores = detectCores())
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
        file_extension <- sub(
            pattern = ".csv", replacement = "",
            x = sub(
                pattern = "metadata_train_classes:",
                replacement = "", x = basename(fmeta_train)
            )
        )

        # Train models
        train_helper <- function(modeli) {
            cat("\nTraining", modeli, "\n")
            file_modeli <- paste0(
                path2save, modeli,
                "_", transformationi,
                "_", file_extension,
                ".RData"
            )
            if (file.exists(file_modeli)) {
                cat("Loading model...\n")
                load(file_modeli)
                return(modeli_fit)
            }
            if (modeli == "svmRadial") {
                modeli_fit <- caret::train(
                    x = datain$train_x, y = datain$train_y,
                    method = "svmRadial",
                    trControl = fit_control,
                    tuneGrid = expand.grid(
                        sigma = kernlab::sigest(
                            datain$train_x,
                            na.action = na.omit, scaled = TRUE
                        ),
                        C = c(0.25, 0.5, 2^c(1:10))
                    ),
                    metric = "ROC",
                    verbose = FALSE
                )
            } else {
                modeli_fit <- caret::train(
                    x = datain$train_x, y = datain$train_y,
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
                    path2save, modeli, "_tuning_",
                    transformationi, "_", file_extension,
                    ".pdf"
                ),
                width = 10, height = 5
            )
            print(
                plot(
                    modeli_fit,
                    ylim = c(
                        floor(
                            min(modeli_fit$results[, "ROC"]) * 10
                        ) / 10, 1
                    ),
                    main = paste(modeli, transformationi, sep = " | ")
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
                path2save, "models_comparison_",
                transformationi,
                "_", file_extension,
                ".RData"
            )
        )
        pdf(
            file = paste0(
                path2save, "models_comparison_", transformationi,
                "_", file_extension,
                ".pdf"
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
}
