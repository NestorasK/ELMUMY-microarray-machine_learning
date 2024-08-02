rm(list = ls())
library(ggplot2)
library(data.table)
library(dplyr)
library(boot)
library(httpgd)
hgd()
hgd_browse()

# Read predictions ####
predictions <- fread(
    file = paste0(
        "results/experiments_caret/multiple_myeloma_progression/",
        "predictions_optimizing_auc/predictions.csv"
    )
)
# Remove because we cannot calculate auc with one class
predictions <- predictions[!dataset %in% c("GSE5900", "GSE14230"), ]
predictions <- predictions[
    !(meta_train == "metadata_train_classes:['MGUS', 'MM']_dataset:['GSE6477', 'GSE2113', 'EMTAB316', 'GSE13591'].csv" &
        dataset == "EMTAB317"),
]
predictions$class[predictions$class == "progressing_MGUS"] <- "MM"

# Path to save
path2save <- paste0(
    "results/experiments_caret/multiple_myeloma_progression/",
    "evaluation_optimizing_auc/"
)

# Collect cv auc ####
path2read_models <- paste0(
    "results/experiments_caret/multiple_myeloma_progression/",
    "training_optimizing_auc"
)

transf_metatrain <- unique(
    predictions[, c("transformation", "meta_train")]
)
auc_cv <- vector(
    mode = "list",
    length = nrow(transf_metatrain)
)
for (i in seq_len(length.out = length(auc_cv))) {
    transf_metatraini <- transf_metatrain[i]
    transformationi <- transf_metatraini[, transformation]
    meta_traini <- transf_metatraini[, meta_train]
    if (transformationi == "qnorm") {
        file_cv_perf <- paste0(
            path2read_models,
            "/models_comparison_",
            transformationi,
            sub(
                pattern = "metadata_train_",
                replacement = "_",
                x = transf_metatraini[, meta_train]
            ),
            collapse = ""
        )
        file_cv_perf <- paste0(
            file_cv_perf, sub(
                pattern = "metadata_train_classes:",
                replacement = "_",
                x = transf_metatraini[, meta_train]
            ), ".RData"
        ) %>% gsub(pattern = ".csv", replacement = "", fixed = TRUE)
    } else {
        file_cv_perf <- paste0(
            path2read_models,
            "/models_comparison_",
            transformationi,
            sub(
                pattern = "metadata_train_classes:",
                replacement = "_", x = transf_metatraini[, meta_train]
            ),
            ".RData"
        ) %>% sub(pattern = ".csv", replacement = "")
    }
    load(file = file_cv_perf)
    auc_dt <- data.table(
        models_comparison$statistics$ROC,
        keep.rownames = TRUE
    )
    auc_dt$transformation <- transformationi
    auc_dt$meta_train <- meta_traini
    auc_cv[[i]] <- auc_dt
}
auc_cv_dt <- do.call(what = "rbind", args = auc_cv)
colnames(auc_cv_dt)[1] <- "method"

# Calculate test auc ####
auc_holdout <- predictions[
    , .(
        auc_test = pROC::auc(
            response = class,
            predictor = MGUS
        )
    ),
    by = c("dataset", "method", "transformation", "meta_train")
]

# Compare ####
auc_all <- merge.data.table(
    x = auc_holdout, y = auc_cv_dt,
    by = c("method", "transformation", "meta_train")
)
auc_all[
    ,
    diff_multiauccv_multiauc_holdout := Mean - auc_test
]
# Platform information
auc_all$platform[
    auc_all$dataset %in% c("GSE2113", "GSE6477", "GSE13591", "GSE14230")
] <- "GLP96"
auc_all$platform[
    auc_all$dataset %in% c("GSE5900", "GSE235356")
] <- "GLP570"
auc_all$platform[
    auc_all$dataset %in% c("EMTAB316")
] <- "A.AFFY.34"
auc_all$platform[
    auc_all$dataset %in% c("EMTAB317")
] <- "A.AFFY.44"


# Performance ####
# All datasets
myfun <- function(str) {
    out <- sub(
        pattern = "metadata_train_classes:['MGUS', 'MM']_dataset:['",
        replacement = "", x = str, fixed = TRUE
    ) %>% sub(pattern = "'].csv", replacement = "", fixed = TRUE)
    out <- paste(
        unlist(
            strsplit(x = out, split = "', '", fixed = TRUE)
        ),
        collapse = " + "
    )
    out <- paste("train:", out)
    return(out)
}
auc_all$meta_train2plot <- sapply(X = auc_all[, meta_train], FUN = myfun)
fwrite(x = auc_all, file = paste0(path2save, "auc_all.csv"))

for (metatri in unique(auc_all[, meta_train2plot])) {
    if (grepl(pattern = "+", x = metatri, fixed = TRUE)) {
        widthi <- 5
        heighti <- 6
    } else {
        widthi <- 5
        heighti <- 8
    }
    plot_jitter_perf <- ggplot(
        data = auc_all[meta_train2plot == metatri, ],
        mapping = aes(
            x = transformation, y = auc_test,
            colour = method
        )
    ) +
        geom_jitter(width = 0.1, height = 0, shape = 1) +
        scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
        ylab("AUC holdout") +
        theme(
            axis.text.x = element_text(angle = 60, hjust = 1),
            legend.position = "bottom"
        ) +
        facet_wrap(
            facets = vars(paste(dataset, platform, sep = " - ")),
            ncol = 2
        ) +
        ggtitle(paste(metatri))
    ggsave(
        filename = paste0(
            path2save, metatri,
            "_plot_jitter_performance_alldata.pdf"
        ),
        plot = plot_jitter_perf, width = widthi, height = heighti
    )
}

for (metatri in unique(auc_all[, meta_train2plot])) {
    plot_point_perf <- ggplot(
        data = melt.data.table(
            data = unique(
                auc_all[
                    meta_train2plot == metatri & dataset == "GSE235356",
                    c(
                        "method", "transformation", "meta_train",
                        "dataset", "platform", "meta_train2plot",
                        "auc_test", "Mean"
                    )
                ]
            ),
            id.vars = c(
                "method", "transformation", "meta_train",
                "dataset", "platform", "meta_train2plot"
            ),
            variable.name = "metric", value.name = "AUC"
        ),
        mapping = aes(
            x = transformation, y = AUC,
            colour = metric
        )
    ) +
        geom_point(shape = 1) +
        scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
        ylab("AUC holdout") +
        theme(
            axis.text.x = element_text(angle = 60, hjust = 1),
            legend.position = "bottom"
        ) +
        facet_wrap(facets = vars(method), ncol = 2) +
        ggtitle(paste0(metatri, " | test: GSE235356"))

    ggsave(
        filename = paste0(
            path2save, metatri,
            "_plot_point_performance_GSE235356.pdf"
        ),
        plot = plot_point_perf, width = 5, height = 8
    )
}
