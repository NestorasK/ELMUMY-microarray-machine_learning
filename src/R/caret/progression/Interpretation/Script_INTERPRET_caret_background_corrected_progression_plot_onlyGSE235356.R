rm(list = ls())
library(data.table)
library(ggplot2)
library(pheatmap)
library(org.Hs.eg.db)
library(annotate)
library(hgu133plus2.db)
library(httpgd)
hgd()
hgd_browse()

path2read <- paste0(
    "results/experiments_caret/multiple_myeloma_progression_onlyGSE235356/",
    "interpretation_optimizing_auc/variable_importance/"
)

importance_all <- fread(
    file = paste0(
        path2read,
        "variable_importance_caret.csv"
    )
)

# Split ratios and not ratios
importance_all_ratios <- importance_all[
    transformation == "ratios" & transformation != "binary_0",
]
importance_all_notratios <- importance_all[
    transformation != "ratios" & transformation != "binary_0",
]
# Add gene information
feature_data <- fread("data/raw/GSE235356/feature_data.csv")
feature_data$Gene.symbol <- sub(
    pattern = "///", replacement = "_",
    x = feature_data$Gene.symbol
)
feature_data$Gene.symbol.short <- substr(
    x = feature_data$Gene.symbol, start = 1, stop = 15
)
feature_data$Gene.symbol.short <- paste0(
    feature_data$Gene.symbol.short,
    ifelse(test = nchar(feature_data$Gene.symbol) > 15,
        yes = "...", no = ""
    )
)

importance_all_notratios <- merge.data.table(
    x = feature_data[, c("rn", "Gene.symbol.short")],
    y = importance_all_notratios
)
importance_all_notratios[, rn2plot := paste0(rn, "(", Gene.symbol.short, ")")]

# Testing - TODO
importance_allj <- copy(importance_all_notratios)

# Mean overall #
importance_allj[
    , myOverall := rowMeans(x = .SD, na.rm = TRUE),
    .SDcols = c("MGUS", "progressing_MGUS", "Overall")
]
importance_allj[, method_transform := paste0(method, "_", transformation)]
data2plot_mean <- unique(
    importance_allj[, c(
        "rn2plot", "method", "transformation", "method_transform",
        "myOverall"
    )]
)

# Calculate ranks ####
# the highest the more important
importance_all_melt <- melt.data.table(
    data = importance_all_notratios[
        , -c("Gene.symbol", "Gene.symbol.short", "rn")
    ],
    id.vars = c("rn2plot", "method", "transformation"),
    variable.name = "perclass", value.name = "importance"
)
importance_all_melt <- na.omit(importance_all_melt)
importance_all_melt[
    , importance_rank_perclass := frank(x = importance),
    by = c("method", "transformation", "perclass")
]
importance_sumrank_dt <- importance_all_melt[
    , .(importance_sumrank = sum(importance_rank_perclass)),
    by = c("rn2plot", "method", "transformation")
]
importance_sumrank_dt[
    , final_rank := frank(
        x = importance_sumrank,
        ties.method = "min"
    ),
    by = c("method", "transformation")
]
importance_sumrank_dt[, method_transform := paste0(method, "_", transformation)]
data2plot_rank <- unique(
    importance_sumrank_dt[, c(
        "rn2plot", "method", "transformation",
        "method_transform", "final_rank"
    )]
)

# Merge both ranks and myOverall
data2plot <- merge.data.table(
    x = data2plot_rank, y = data2plot_mean,
    by = c("rn2plot", "method", "transformation", "method_transform")
)

# Relation between ranks and myOverall
plot(
    data2plot[
        sample(x = 1:nrow(data2plot), size = 10000),
        c("final_rank", "myOverall")
    ],
    main = paste(
        "spearman:",
        round(
            cor(data2plot$final_rank, data2plot$myOverall, method = "sp"),
            digits = 2
        )
    )
)

# - Per method score ####
# scorei <- c("final_rank", "myOverall")[1]
# methodi <- unique(data2plot[, method])[1]

# topnumgenes <- length(unique(data2plot[, rn]))
topnumgenes <- "all"
topnumgenes <- 50
for (scorei in c("final_rank", "myOverall")) {
    cat("Plotting:", scorei, "\n")
    for (methodi in unique(data2plot[, method])) {
        cat(" -", methodi, "\n")
        importance_dcastj <- dcast.data.table(
            data = data2plot[method == methodi, ],
            formula = rn2plot ~ transformation,
            value.var = scorei
        )
        # Keep X most important probes
        importance_dcastj$sum_importance <- rowSums(
            as.matrix(importance_dcastj[, -"rn2plot"])
        )
        importance_dcastj <- importance_dcastj[order(-sum_importance), ]
        if (topnumgenes == "all") {
            if (scorei == "myOverall") {
                topnumgenesi <- sum(importance_dcastj$sum_importance > 0)
            } else if (scorei == "final_rank") {
                topnumgenesi <- sum(importance_dcastj$sum_importance > 4)
            }
        } else {
            topnumgenesi <- topnumgenes
        }
        importance_dcast_mat <- as.matrix(
            importance_dcastj[1:topnumgenesi, -c("rn2plot", "sum_importance")]
        )
        rownames(importance_dcast_mat) <- importance_dcastj[
            1:topnumgenesi, rn2plot
        ]
        my_heatmap <- pheatmap(
            mat = importance_dcast_mat,
            color = colorRampPalette(c("#fffec3", "darkred"))(10000),
            show_rownames = TRUE,
            fontsize_col = 10,
            fontsize_row = 7,
            main = paste0(
                "Top ", topnumgenesi, " important features | ",
                methodi
            )
        )
        png(
            filename = paste0(
                path2read,
                "heatmap_variable_importance_",
                methodi, "_", scorei, "_topnumbergenes_", topnumgenesi,
                "_onlyGSE235356.png"
            ), width = 1300, height = 2400,
            res = 300
        )
        grid::grid.newpage()
        grid::grid.draw(my_heatmap$gtable)
        dev.off()
    }
    # All together
    cat(" - all together\n")
    importance_dcast <- dcast.data.table(
        data = data2plot[method != "svmRadial", ],
        formula = rn2plot ~ method_transform,
        value.var = scorei
    )
    # Keep X most important probes
    importance_dcast$sum_importance <- rowSums(
        as.matrix(importance_dcast[, -"rn2plot"])
    )
    importance_dcast <- importance_dcast[order(-sum_importance), ]
    if (topnumgenes == "all") {
        if (scorei == "myOverall") {
            topnumgenesi <- sum(importance_dcastj$sum_importance > 0)
        } else if (scorei == "final_rank") {
            topnumgenesi <- sum(importance_dcastj$sum_importance > 4)
        }
    } else {
        topnumgenesi <- topnumgenes
    }
    importance_dcast_mat <- as.matrix(
        importance_dcast[1:topnumgenesi, -c("rn2plot", "sum_importance")]
    )
    rownames(importance_dcast_mat) <- importance_dcast[1:topnumgenesi, rn2plot]

    my_heatmap_all <- pheatmap(
        mat = importance_dcast_mat,
        color = colorRampPalette(c("#fffec3", "darkred"))(30000),
        show_rownames = TRUE,
        fontsize_col = 10,
        fontsize_row = 7,
        main = paste0(
            "Top ", topnumgenes, " important features | all methods"
        )
    )
    png(
        filename = paste0(
            path2read,
            "heatmap_variable_importance_allmethods_",
            scorei, "_topnumbergenes_", topnumgenesi,
            "_onlyGSE235356.png"
        ), width = 1600, height = 2400,
        res = 300
    )
    grid::grid.newpage()
    grid::grid.draw(my_heatmap_all$gtable)
    dev.off()
}
