rm(list = ls())
library(data.table)
library(ggplot2)
library(pheatmap)
library(httpgd)
hgd()
hgd_browse()

path2read <- paste0(
    "results/experiments_caret/multiple_myeloma_stage/",
    "interpretation_optimizing_multiclass_auc/"
)

importance_all <- fread(
    file = paste0(
        path2read,
        "variable_importance_caret.csv"
    )
)

importance_all[
    , myOverall := rowMeans(x = .SD, na.rm = TRUE),
    .SDcols = c("MGUS", "MM", "Normal", "Overall")
]
numfeatures <- importance_all[
    transformation != "binary_0",
    .(nfeatures = sum(myOverall > 0)),
    by = c("method", "transformation")
]
mean_nfeatures <- numfeatures[, .(mean_methods = mean(nfeatures)), by = method]
numfeatures$method <- factor(
    numfeatures$method,
    levels = mean_nfeatures[
        order(mean_methods, na.last = TRUE, decreasing = FALSE), method
    ]
)

dot_plot_nimpfeat <- ggplot(
    data = numfeatures,
    mapping = aes(x = method, y = nfeatures, colour = transformation)
) +
    geom_jitter(width = 0.1, height = 0, shape = 1) +
    xlab("") +
    ylab("Number of features") +
    ggtitle(
        label = "Number of features in the models"
    ) +
    coord_flip()

ggsave(
    filename = paste0(
        "results/experiments_caret/multiple_myeloma_stage/",
        "interpretation_optimizing_multiclass_auc/",
        "dotplot_num_imp_features_myOverall.pdf"
    ), plot = dot_plot_nimpfeat, width = 5, height = 2.5
)
