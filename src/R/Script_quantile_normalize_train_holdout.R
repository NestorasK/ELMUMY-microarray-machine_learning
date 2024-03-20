rm(list = ls())
library(data.table)
library(preprocessCore)
library(ggplot2)

# Input ####
path2read <- "data/processed_gpl96_gpl570_affy44_platform/"
file_expr <- paste0(
    path2read,
    "background_corrected_expression_rma.csv"
)
fsmeta_train <- list.files(
    path = path2read, pattern = "metadata_train_",
    full.names = TRUE, recursive = TRUE
)
fsmeta_test <- sub(
    pattern = "metadata_train_", replacement = "metadata_holdout_",
    x = fsmeta_train
)

# Calculations ####
expr <- fread(file_expr)
ifile <- 1
for (ifile in seq_len(length.out = length(fsmeta_train))) {
    cat("\nPerforming quantile normalization on:", fsmeta_train[ifile], "\n")

    meta_train <- fread(fsmeta_train[ifile])
    meta_holdout <- fread(fsmeta_test[ifile])

    path2save <- paste(
        strsplit(
            x = fsmeta_train[ifile], split = "/", fixed = TRUE
        )[[1]][c(1, 2)],
        collapse = "/"
    )

    source("src/R/fetch_data.R")
    datain <- fetch_data(
        expres = expr,
        meta_train = meta_train, meta_holdout = meta_holdout
    )
    # train
    train_x_qnorm <- t(normalize.quantiles(t(datain$train_x)))
    rownames(train_x_qnorm) <- rownames(datain$train_x)
    colnames(train_x_qnorm) <- colnames(datain$train_x)

    # holdout
    holdout_x_qnorm <- t(
        normalize.quantiles.use.target(
            x = t(datain$holdout_x),
            target = train_x_qnorm[1, ]
        )
    )
    rownames(holdout_x_qnorm) <- rownames(datain$holdout_x)
    colnames(holdout_x_qnorm) <- colnames(datain$holdout_x)

    # merge data
    merge_x_qnorm <- merge.data.table(
        x = data.table(t(train_x_qnorm), keep.rownames = TRUE),
        y = data.table(t(holdout_x_qnorm), keep.rownames = TRUE),
        by = "rn"
    )
    filename_qnorm <- sub(
        pattern = "metadata_train_",
        replacement = "background_corrected_expression_qnorm_",
        x = fsmeta_train[ifile]
    )
    cat("Write qnorm file:", filename_qnorm, "\n")
    fwrite(x = merge_x_qnorm, file = filename_qnorm)

    # diagnostic plots
    filename_diagnostic_plots <- sub(
        pattern = ".csv", replacement = ".pdf",
        x = sub(
            pattern = "metadata_train_",
            replacement = "diagnostic_plot_quantile_norm_",
            x = fsmeta_train[ifile]
        )
    )
    cat("Write diagnostics plot file:", filename_diagnostic_plots, "\n")
    pdf(
        file = filename_diagnostic_plots,
        width = 5, height = 10
    )
    par(mfrow = c(2, 1))
    df_background <- as.data.frame(
        t(rbind(datain$train_x[1:10, ], datain$holdout_x[1:10, ]))
    )
    colnames(df_background) <- NULL
    boxplot(
        df_background,
        col = rep(c("grey", "white"), each = 10),
        las = 1,
        main = paste("background corrected\n", basename(path2save)),
        ylim = ceiling(range(df_background)) - c(5, 0)
    )
    legend(
        x = "bottomright",
        legend = c("Training set", "Holdout"),
        fill = c("grey", "white")
    )
    df_qnorm <- t(rbind(train_x_qnorm[1:10, ], holdout_x_qnorm[1:10, ]))
    colnames(df_qnorm) <- NULL
    boxplot(df_qnorm,
        col = rep(c("grey", "white"), each = 10), las = 1,
        main = paste(
            "background corrected + quantile norm\n",
            basename(path2save)
        ),
        ylim = ceiling(range(df_qnorm)) - c(5, 0)
    )
    legend(
        x = "bottomright",
        legend = c("Training set", "Holdout"),
        fill = c("grey", "white")
    )
    dev.off()
}
