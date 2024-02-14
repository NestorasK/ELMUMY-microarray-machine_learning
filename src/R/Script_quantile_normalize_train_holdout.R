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
    path = "data", pattern = "metadata_train.csv",
    full.names = TRUE, recursive = TRUE
)
fsmeta_test <- sub(
    pattern = "metadata_train.csv", replacement = "metadata_holdout.csv",
    x = fsmeta_train
)

# Calculations ####
expr <- fread(file_expr)

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
    fwrite(
        x = data.table(train_x_qnorm, keep.rownames = TRUE),
        file = paste0(path2save, "/qnorm_expression_train.csv")
    )

    # holdout
    holdout_x_qnorm <- t(
        normalize.quantiles.use.target(
            x = t(datain$holdout_x),
            target = train_x_qnorm[1, ]
        )
    )
    rownames(holdout_x_qnorm) <- rownames(datain$holdout_x)
    colnames(holdout_x_qnorm) <- colnames(datain$holdout_x)
    fwrite(
        x = data.table(holdout_x_qnorm, keep.rownames = TRUE),
        file = paste0(path2save, "/qnorm_expression_holdout.csv")
    )

    # diagnostic plots
    pdf(
        file = paste0(path2save, "/diagnostic_plot_quantile_norm.pdf"),
        width = 10, height = 20
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
