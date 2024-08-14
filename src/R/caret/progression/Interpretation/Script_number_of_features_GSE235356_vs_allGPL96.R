rm(list = ls())
library(data.table)
library(ggplot2)
library(pheatmap)
library(VennDiagram)
library(httpgd)
hgd()
hgd_browse()


path2save <- paste0(
    "results/experiments_caret/multiple_myeloma_progression/",
    "interpretation_optimizing_auc/variable_importance/"
)

# Variable importance - ALL GPL96
path2read <- paste0(
    "results/experiments_caret/multiple_myeloma_progression/",
    "interpretation_optimizing_auc/variable_importance/"
)
importance_all_gpl96 <- fread(
    file = paste0(
        path2read,
        "variable_importance_caret.csv"
    )
)
importance_all_gpl96$meta_train <- "all GLP96"
importance_all_gpl96[
    , myOverall := rowMeans(x = .SD, na.rm = TRUE),
    .SDcols = c("MGUS", "MM", "Overall")
]


# Variable importance - onlyGSE235356
path2read <- paste0(
    "results/experiments_caret/multiple_myeloma_progression_onlyGSE235356/",
    "interpretation_optimizing_auc/variable_importance/"
)
importance_onlyGSE235356 <- fread(
    file = paste0(
        path2read,
        "variable_importance_caret.csv"
    )
)
colnames(importance_onlyGSE235356)
importance_onlyGSE235356$meta_train <- "GSE235356"
importance_onlyGSE235356[
    , myOverall := rowMeans(x = .SD, na.rm = TRUE),
    .SDcols = c("MGUS", "progressing_MGUS", "Overall")
]

# Merge
importance_all <- rbindlist(
    l = list(
        importance_all_gpl96[
            , c("rn", "method", "transformation", "meta_train", "myOverall")
        ],
        importance_onlyGSE235356[
            , c("rn", "method", "transformation", "meta_train", "myOverall")
        ]
    )
)
numfeatures <- importance_all[
    transformation != "binary_0",
    .(nfeatures = sum(myOverall > 0)),
    by = c("method", "transformation", "meta_train")
]

tapply(
    X = numfeatures$nfeatures,
    INDEX = paste(numfeatures$meta_train, numfeatures$method),
    FUN = summary
)

set.seed(123)
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
    # coord_flip() +
    theme(
        axis.text.x = element_text(angle = 60, hjust = 1),
        legend.position = "bottom"
    ) +
    facet_wrap(facets = vars(meta_train))
ggsave(
    filename = paste0(
        path2save,
        "dotplot_num_imp_features_myOverall_onlyGSE_vs_allGPL96.pdf"
    ),
    plot = dot_plot_nimpfeat, width = 6, height = 4
)

# Number of common features
number_of_common_probes <- unique(
    importance_all[
        (transformation != "binary_0") &
            (method %in% c("glmnet", "rf", "gbm")) &
            (myOverall > 0), .(number_of_methods = length(unique(meta_train))),
        by = c("transformation", "rn", "method")
    ]
)
number_of_common_probes[, method_transf := paste0(method, "+", transformation)]

number_of_methods <- rbindlist(
    l = lapply(
        X = unique(number_of_common_probes$method_transf),
        FUN = function(trj, mati = number_of_common_probes) {
            dt <- data.table(
                table(mati[method_transf == trj, number_of_methods])
            )
            dt$method_transf <- trj
            return(dt)
        }
    )
)
colnames(number_of_methods)[c(1, 2)] <- c(
    "number_of_trainingsets", "number_of_probes"
)
number_of_methods[
    , c("method", "transformation") := tstrsplit(
        x = method_transf,
        split = "+", keep = c(1, 2), fixed = TRUE
    )
]

number_of_methods <- number_of_methods[number_of_trainingsets == 2, ]
tapply(
    X = number_of_methods$number_of_probes,
    INDEX = paste(number_of_methods$method),
    FUN = summary
)

set.seed(42)
dotplot_commonprobes_methods <- ggplot(
    data = number_of_methods[number_of_trainingsets == 2, ],
    mapping = aes(
        x = method,
        y = number_of_probes, colour = transformation
    )
) +
    geom_jitter(width = 0.1, height = 0, shape = 1) +
    ylab("Number of common probes") +
    xlab("") +
    ggtitle(
        label = "Number of common probes across training sets",
        subtitle = "all GPL96 | GSE235356"
    )
ggsave(
    filename = paste0(
        path2save,
        "Number of common probes across onlyGSE_vs_allGPL96.pdf"
    ),
    plot = dotplot_commonprobes_methods,
    width = 6, height = 3.5
)

# Venn diagrams ####
for (trj in unique(importance_all[, transformation])) {
    dtj <- dcast.data.table(
        data = importance_all[(transformation == trj) & (myOverall != 0)],
        formula = rn + method ~ meta_train,
        value.var = "myOverall",
        fun.aggregate = length
    )
    for (mi in unique(dtj[, method])) {
        dtj_mi <- dtj[method == mi, ]
        vennplot <- venn.diagram(
            x = list(
                GSE235356 = which(dtj_mi$GSE235356 > 0),
                GLP96 = which(dtj_mi$all_GLP96 > 0)
            ),
            filename = NULL,
            main = paste0(
                "Number of common probes across training sets\n",
                "Method:", mi, " | Transformations:", trj
            )
        )
        pdf(
            file = paste0(
                path2save,
                "venn_diagram_common_features_across_trainingsets_",
                trj, "_transformation_",
                mi, "_method",
                ".pdf"
            ),
            width = 4.5, height = 4
        )
        grid.draw(vennplot)
        dev.off()
    }
}
