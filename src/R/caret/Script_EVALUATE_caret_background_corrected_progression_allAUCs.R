rm(list = ls())
library(ggplot2)
library(data.table)
library(dplyr)
library(httpgd)
hgd()
hgd_browse()

auc_notGSE235356 <- fread(
    paste0(
        "results/experiments_caret/multiple_myeloma_progression/",
        "evaluation_optimizing_auc/auc_all.csv"
    )
)
auc_notGSE235356$ncv_fold <- 1
auc_GSE235356 <- fread(
    paste0(
        "results/experiments_caret/multiple_myeloma_progression_onlyGSE235356/",
        "evaluation_optimizing_auc/auc_all.csv"
    )
)
auc_GSE235356$dataset <- "GSE235356"
auc_GSE235356$platform <- "GLP570"
auc_GSE235356$meta_train2plot <- "train: GSE235356"

auc_all <- rbindlist(
    l = list(
        auc_notGSE235356[
            transformation != "binary_0" & dataset == "GSE235356",
            c(
                "method", "transformation", "dataset",
                "auc_test", "Mean",
                "meta_train2plot"
            )
        ],
        auc_GSE235356[
            , c(
                "method", "transformation", "dataset",
                "auc_test", "Mean",
                "meta_train2plot"
            )
        ]
    ), use.names = TRUE
)

colnames(auc_all)[6] <- "training"
auc_all$training <- sub(
    pattern = "train: ", replacement = "",
    x = auc_all$training
)
plot_auc <- ggplot(
    data = auc_all[training == "GSE235356"],
    mapping = aes(x = transformation, y = auc_test)
) +
    geom_boxplot() +
    geom_point(
        data = auc_all[training != "GSE235356"],
        mapping = aes(
            x = transformation, y = auc_test,
            colour = training
        )
    ) +
    ylab("AUC holdout") +
    xlab("") +
    # theme_bw() +
    theme(
        axis.text.x = element_text(angle = 60, hjust = 1),
        legend.position = "bottom"
    ) +
    facet_wrap(facets = vars(method), ncol = 2)
ggsave(
    filename = paste0(
        "results/experiments_caret/multiple_myeloma_progression_onlyGSE235356/",
        "evaluation_optimizing_auc/plot_generalization_performance_alldata.pdf"
    ), width = 6, height = 8
)


# Plot using GSE235356 cross validation error as background
auc_GSE235356[, auc_test_mean := mean(auc_test),
    by = c("method", "transformation")
]

auc_GSE235356 <- unique(
    x = auc_GSE235356[
        , c(
            "method", "transformation", "dataset", "meta_train2plot", "Mean",
            "auc_test_mean"
        )
    ]
)
colnames(auc_GSE235356)[6] <- "auc_test"

auc_all <- rbindlist(
    l = list(
        auc_notGSE235356[
            transformation != "binary_0" & dataset == "GSE235356",
            c(
                "method", "transformation", "dataset",
                "auc_test", "Mean",
                "meta_train2plot"
            )
        ],
        auc_GSE235356[
            , c(
                "method", "transformation", "dataset",
                "auc_test", "Mean",
                "meta_train2plot"
            )
        ]
    ), use.names = TRUE
)

colnames(auc_all)[6] <- "training"
auc_all$training <- sub(
    pattern = "train: ", replacement = "",
    x = auc_all$training
)

plot_auccv <- ggplot(
    data = auc_all[training == "GSE235356"],
    mapping = aes(x = transformation, y = Mean)
) +
    geom_boxplot() +
    geom_jitter(
        data = unique(auc_all[, -c("Mean")]),
        mapping = aes(
            x = transformation, y = auc_test,
            colour = training
        ), width = 0.1, height = 0
    ) +
    ylab("AUC") +
    xlab("") +
    ggtitle("Test set: GSE235356") +
    # theme_bw() +
    theme(
        axis.text.x = element_text(angle = 60, hjust = 1),
        legend.position = "bottom"
    ) +
    facet_wrap(facets = vars(method), ncol = 2)
ggsave(
    filename = paste0(
        "results/experiments_caret/multiple_myeloma_progression_onlyGSE235356/",
        "evaluation_optimizing_auc/plot_generalization_performance_alldata_cverror.pdf"
    ), plot = plot_auccv,
    width = 6, height = 9
)
