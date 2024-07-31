rm(list = ls())
library(data.table)
library(ggplot2)
library(pheatmap)
library(VennDiagram)
library(httpgd)
hgd()
hgd_browse()

path2read <- paste0(
    "results/experiments_caret/multiple_myeloma_stage/",
    "interpretation_optimizing_multiclass_auc/variable_importance/"
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
        "interpretation_optimizing_multiclass_auc/variable_importance/",
        "dotplot_num_imp_features_myOverall.pdf"
    ), plot = dot_plot_nimpfeat, width = 5, height = 2.5
)

# Number of common features
number_of_common_probes <- unique(
    importance_all[
        (transformation != "binary_0") &
            (method %in% c("glmnet", "rf", "gbm")) &
            (myOverall > 0), .(number_of_methods = length(unique(method))),
        by = c("transformation", "rn")
    ]
)
number_of_methods <- rbindlist(
    l = lapply(
        X = unique(number_of_common_probes$transformation),
        FUN = function(trj, mati = number_of_common_probes) {
            dt <- data.table(table(mati[transformation == trj, number_of_methods]))
            dt$transformation <- trj
            colnames(dt)[c(1, 2)] <- c("number_of_methods", "number_of_probes")
            return(dt)
        }
    )
)
set.seed(123)
dotplot_commonprobes_methods <- ggplot(
    data = number_of_methods,
    mapping = aes(
        x = number_of_methods,
        y = number_of_probes, colour = transformation
    )
) +
    geom_jitter(width = 0.2, height = 0) +
    ylab("Number of probes") +
    xlab("Number of methods") +
    ggtitle(
        label = "Number of common probes across ML methods",
        subtitle = "Methods: glmnet, gbm, rf"
    )
ggsave(
    filename = paste0(
        "results/experiments_caret/multiple_myeloma_stage/",
        "interpretation_optimizing_multiclass_auc/variable_importance/",
        "Number of common probes across ML methods.pdf"
    ),
    plot = dotplot_commonprobes_methods,
    width = 6, height = 3.5
)

# Number of common probes across transformations
probes_info <- unique(
    importance_all[myOverall > 0, c("rn", "method", "transformation")]
)
probes_ratios <- rbindlist(
    l = lapply(
        X = unique(probes_info[, method]),
        FUN = function(mi) {
            data.table(
                rn = unlist(
                    strsplit(
                        x = probes_info[
                            transformation == "ratios" &
                                method == mi, rn
                        ],
                        split = "///",
                        fixed = TRUE
                    )
                ),
                method = mi
            )
        }
    )
)
probes_ratios$transformation <- "ratios"
probes_info_all <- rbindlist(
    l = list(
        unique(probes_info[transformation != "ratios"]),
        unique(probes_ratios)
    )
)
number_of_common_probes <- unique(
    probes_info_all[
        (transformation != c("binary_0")) &
            (method %in% c("glmnet", "rf", "gbm")), .(
            number_of_transformation = length(unique(transformation))
        ),
        by = c("method", "rn")
    ]
)
number_of_transformation_all <- rbindlist(
    l = lapply(
        X = unique(number_of_common_probes$method),
        FUN = function(trj, mati = number_of_common_probes) {
            dt <- data.table(table(mati[
                method == trj,
                number_of_transformation
            ]))
            dt$method <- trj
            colnames(dt)[c(1, 2)] <- c(
                "number_of_transformation", "number_of_probes"
            )
            return(dt)
        }
    )
)
set.seed(42)
dotplot_commonprobes_transf <- ggplot(
    data = number_of_transformation_all,
    mapping = aes(
        x = number_of_transformation,
        y = number_of_probes, colour = method
    )
) +
    geom_jitter(width = 0.2, height = 0) +
    ylab("Number of probes") +
    xlab("Number of transformations") +
    ggtitle(
        label = "Number of common probes across transformations",
        subtitle = "Transformations: rma, binary_0.5, qnorm, ranking, ratios"
    )
ggsave(
    filename = paste0(
        "results/experiments_caret/multiple_myeloma_stage/",
        "interpretation_optimizing_multiclass_auc/variable_importance/",
        "Number of common probes across transformations.pdf"
    ),
    plot = dotplot_commonprobes_transf,
    width = 6, height = 3.5
)
# Plot together
set.seed(123)
ggsave(
    filename = paste0(
        "results/experiments_caret/multiple_myeloma_stage/",
        "interpretation_optimizing_multiclass_auc/variable_importance/",
        "Number of common probes across.pdf"
    ),
    plot = cowplot::plot_grid(
        plotlist = list(
            dotplot_commonprobes_methods,
            dotplot_commonprobes_transf
        ),
        nrow = 2, labels = "AUTO"
    ),
    width = 6, height = 6
)

# Venn diagrams
for (trj in unique(importance_all[, transformation])) {
    dtj <- dcast.data.table(
        data = importance_all[(transformation == trj) & (myOverall != 0)],
        formula = rn ~ method,
        value.var = "myOverall",
        fun.aggregate = length
    )
    vennplot <- venn.diagram(
        x = list(
            gbm = which(dtj$gbm > 0),
            glmnet = which(dtj$glmnet > 0),
            rf = which(dtj$rf > 0)
        ),
        filename = NULL,
        main = paste0("Number of common probes across ML methods | ", trj)
    )
    pdf(
        file = paste0(
            "results/experiments_caret/multiple_myeloma_stage/",
            "interpretation_optimizing_multiclass_auc/variable_importance/",
            "venn_diagram_common_features_across_methods_",
            trj, "_transformation.pdf"
        ),
        width = 4.5, height = 4
    )
    grid.draw(vennplot)
    dev.off()
}

for (trj in unique(probes_info_all[, method])) {
    dtj <- dcast.data.table(
        data = probes_info_all[(method == trj) & (transformation != "binary_0")],
        formula = rn ~ transformation,
        fun.aggregate = length
    )
    vennplot <- venn.diagram(
        x = list(
            binary_0.5 = which(dtj$binary_0.5 > 0),
            qnorm = which(dtj$qnorm > 0),
            ranking = which(dtj$ranking > 0),
            rma = which(dtj$rma > 0),
            ratios = which(dtj$ratios > 0)
        ),
        filename = NULL,
        main = paste0("Number of common probes across transformations | ", trj)
    )
    pdf(
        file = paste0(
            "results/experiments_caret/multiple_myeloma_stage/",
            "interpretation_optimizing_multiclass_auc/variable_importance/",
            "venn_diagram_common_features_across_transforation_",
            trj, "_method.pdf"
        ),
        width = 5, height = 5
    )
    grid.draw(vennplot)
    dev.off()
}
