rm(list = ls())
library(ggplot2)
library(data.table)
library(cowplot)
library(boot)
library(httpgd)
hgd()
hgd_browse()

# Read predictions ####
predictions <- fread(
    file = paste0(
        "results/experiments_caret/multiple_myeloma_stage/",
        "predictions_optimizing_multiclass_auc/predictions.csv"
    )
)

# Path to save
path2save <- paste0(
    "results/experiments_caret/multiple_myeloma_stage/",
    "evaluation_optimizing_multiclass_auc/"
)

# Collect cv multiclass_auc ####
path2read_models <- paste0(
    "results/experiments_caret/multiple_myeloma_stage/",
    "training_optimizing_multiclass_auc/"
)
predictions$transformation[
    grepl(pattern = "qnorm", x = predictions$transformation)
] <-
    "qnorm_classes:['Normal', 'MGUS', 'MM']_dataset:['GSE6477']"
multiclass_auc_cv <- vector(
    mode = "list",
    length = length(unique(predictions$transformation))
)
for (i in seq_len(length.out = length(unique(predictions$transformation)))) {
    transformationi <- unique(predictions$transformation)[i]
    load(
        file = paste0(
            path2read_models, "models_comparison_",
            transformationi, ".RData"
        )
    )
    multiauc_df <- data.table(models_comparison$statistics$multiclass_auc,
        keep.rownames = TRUE
    )
    multiauc_df$transformation <- transformationi
    multiclass_auc_cv[[i]] <- multiauc_df
}
multiclass_auc_cv_dt <- do.call(what = "rbind", args = multiclass_auc_cv)
colnames(multiclass_auc_cv_dt)[1] <- "method"

# Calculate test multiclass_auc ####
predictions$class[predictions$class == "progressing_MGUS"] <- "MM"
myfun <- function(classin, MGUS, MM, Normal) {
    pROC::auc(
        pROC::multiclass.roc(
            response = classin,
            predictor = cbind(MGUS, MM, Normal)
        )
    )
}

multiclass_auc_test <- predictions[
    , .(
        multiclass_auc_test = myfun(
            classin = class, MGUS = MGUS, MM = MM,
            Normal = Normal
        )
    ),
    by = c("dataset", "method", "transformation")
]

# Compare ####
multiclass_auc_all <- merge.data.table(
    x = multiclass_auc_test, y = multiclass_auc_cv_dt,
    by = c("method", "transformation")
)
multiclass_auc_all[
    ,
    diff_multiauccv_multiauc_test := Mean - multiclass_auc_test
]
multiclass_auc_all$transformation[
    multiclass_auc_all$transformation ==
        "qnorm_classes:['Normal', 'MGUS', 'MM']_dataset:['GSE6477']"
] <- "qnorm"

# Performance ####
# All datasets
plot_jitter_perf <- ggplot(
    data = multiclass_auc_all[dataset != "GSE235356", ],
    mapping = aes(
        x = transformation, y = multiclass_auc_test,
        colour = method
    )
) +
    geom_jitter(width = 0.1, shape = 1) +
    scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
    ylab("Multiclass AUC (test)") +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom"
    ) +
    facet_wrap(facets = vars(dataset), ncol = 2)
ggsave(
    filename = paste0(path2save, "plot_jitter_performance_alldata.pdf"),
    plot = plot_jitter_perf, width = 6, height = 8
)

plot_perf_box <- list()
counter <- 1
plot_perf_box[[counter]] <- ggplot(
    data = multiclass_auc_all[
        transformation != "binary_0" & dataset != "GSE235356",
    ],
    mapping = aes(x = dataset, y = multiclass_auc_test)
) +
    geom_boxplot() +
    geom_jitter(width = 0.1, shape = 1, colour = "#666464", size = 0.7) +
    scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
    ylab("Multiclass AUC test") +
    ggtitle("All data") +
    coord_flip()


# Keep datasets from GPL96 platform
counter <- counter + 1
plot_perf_box[[counter]] <- ggplot(
    data = multiclass_auc_all[
        (transformation != "binary_0") &
            (dataset %in% c("EMTAB316", "GSE2113", "GSE13591")),
    ],
    mapping = aes(x = transformation, y = multiclass_auc_test)
) +
    geom_boxplot() +
    geom_jitter(
        shape = 1,
        width = 0.2,
        size = 0.7,
        colour = "#666464",
    ) +
    scale_y_continuous(breaks = seq(0.5, 1, 0.1), limits = c(0.5, 1)) +
    ylab("Multiclass AUC test") +
    ggtitle("GPL96 platform datasets") +
    coord_flip()

counter <- counter + 1
plot_perf_box[[counter]] <- ggplot(
    data = multiclass_auc_all[
        (transformation != "binary_0") &
            (dataset %in% c("EMTAB316", "GSE2113", "GSE13591")),
    ],
    mapping = aes(x = method, y = multiclass_auc_test)
) +
    geom_boxplot() +
    geom_jitter(
        shape = 1,
        width = 0.2,
        size = 0.7,
        colour = "#666464"
    ) +
    scale_y_continuous(breaks = seq(0.5, 1, 0.1), limits = c(0.5, 1)) +
    ylab("Multiclass AUC test") +
    ggtitle("GPL96 platform datasets") +
    coord_flip()

ggsave(
    filename = paste0(
        path2save, "plot_allboxplots_performance.pdf"
    ),
    plot = plot_grid(plotlist = plot_perf_box, labels = "AUTO"),
    width = 8, height = 5
)

# Generalization ####
plot_gen_alldata <- ggplot(
    data = multiclass_auc_all[dataset != "GSE235356", ],
    mapping = aes(
        x = transformation, y = diff_multiauccv_multiauc_holdout,
        colour = method
    )
) +
    geom_jitter(shape = 1, width = 0.1) +
    scale_y_continuous(
        breaks = seq(-0.4, 0.9, 0.2),
        limits = c(-0.4, 0.9)
    ) +
    ylab("CV Multiclass AUC - Test Multiclass AUC") +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom"
    ) +
    facet_wrap(facets = vars(dataset), ncol = 2)
ggsave(
    filename = paste0(
        path2save, "plot_jitter_generalization_alldata.pdf"
    ),
    plot = plot_gen_alldata, width = 6, height = 8
)

plot_allbox_gen <- list()
counter <- 1
plot_allbox_gen[[counter]] <- ggplot(
    data = multiclass_auc_all[
        transformation != "binary_0" & dataset != "GSE235356",
    ],
    mapping = aes(
        x = dataset, y = diff_multiauccv_multiauc_holdout
    )
) +
    geom_boxplot() +
    geom_jitter(width = 0.1, shape = 1, colour = "#666464", size = 0.7) +
    scale_y_continuous(breaks = seq(0, 1, 0.2)) +
    ylab("CV Multiclass AUC - Test Multiclass AUC") +
    ggtitle("All data") +
    coord_flip()
# ggsave(
#     filename = paste0(
#         path2save, "plot_boxplot_generalization_alldata.pdf"
#     ),
#     plot = plot_gen_alldata, width = 5, height = 4
# )

counter <- counter + 1
plot_allbox_gen[[counter]] <- ggplot(
    data = multiclass_auc_all[
        (transformation != "binary_0") &
            (dataset %in% c("EMTAB316", "GSE2113", "GSE13591")),
    ],
    mapping = aes(x = transformation, y = diff_multiauccv_multiauc_holdout)
) +
    geom_boxplot() +
    geom_jitter(
        # data = multiclass_auc_all[
        #     (transformation != "binary_0") &
        #         (dataset %in% c("EMTAB316", "GSE2113", "GSE13591")),
        # ],
        # mapping = aes(
        #     x = transformation, y = diff_multiauccv_multiauc_holdout,
        #     colour = method
        # ),
        shape = 1,
        width = 0.2,
        size = 0.7,
        colour = "#666464"
    ) +
    ylab("CV Multiclass AUC - Test Multiclass AUC") +
    coord_flip() +
    ggtitle("GPL96 platform datasets")
# ggsave(
#     filename = paste0(
#         path2save, "plot_boxplot_generalization_gpl96~transformation.pdf"
#     ),
#     plot = plot_generalization_transf, width = 5, height = 4
# )

counter <- counter + 1
plot_allbox_gen[[counter]] <- ggplot(
    data = multiclass_auc_all[
        (transformation != "binary_0") &
            (dataset %in% c("EMTAB316", "GSE2113", "GSE13591")),
    ],
    mapping = aes(x = method, y = diff_multiauccv_multiauc_holdout)
) +
    geom_boxplot() +
    geom_jitter(
        # data = multiclass_auc_all[
        #     (transformation != "binary_0") &
        #         (dataset %in% c("EMTAB316", "GSE2113", "GSE13591")),
        # ],
        # mapping = aes(
        #     x = method, y = diff_multiauccv_multiauc_holdout,
        #     colour = transformation
        # ),
        shape = 1,
        width = 0.2, size = 0.7,
        colour = "#666464"
    ) +
    ylab("CV Multiclass AUC - Test Multiclass AUC") +
    ggtitle("GPL96 platform datasets") +
    coord_flip()
ggsave(
    filename = paste0(
        path2save, "plot_allboxplot_generalization.pdf"
    ),
    plot = plot_grid(plotlist = plot_allbox_gen, labels = "AUTO"),
    width = 8, height = 5
)
