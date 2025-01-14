rm(list = ls())
library(data.table)
library(ggplot2)
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)
library(annotate)
library(hgu133plus2.db)
library(httpgd)
hgd()
hgd_browse()

path2save <- paste0(
    "results/experiments_caret/multiple_myeloma_progression/",
    "interpretation_optimizing_auc/enrichment_analysis/"
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
importance_all_gpl96$meta_train <- "allGLP96"
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
importance_onlyGSE235356$meta_train <- "GSE"
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
importance_all <- importance_all[order(-myOverall), ]
importance_all <- importance_all[
    transformation != "binary_0" & myOverall > 0 & method != "svmLinear2"
]
method_train <- unique(importance_all[, c("method", "meta_train")])

probe_lists <- vector(
    mode = "list",
    length = nrow(method_train)
)

for (j in seq_len(nrow(method_train))) {
    probe_lists[[j]] <- unique(
        x = unlist(
            strsplit(
                x = importance_all[
                    method == method_train[j, method] &
                        meta_train == method_train[j, meta_train], rn
                ],
                split = "///", fixed = TRUE
            )
        )
    )

    # write file
    dt <- data.table(
        unique(
            x = unlist(
                strsplit(
                    x = importance_all[
                        method == method_train[j, method] &
                            meta_train == method_train[j, meta_train], rn
                    ],
                    split = "///", fixed = TRUE
                )
            )
        )
    )
    colnames(dt) <- paste0(
        ">", method_train[j, method],
        "_", method_train[j, meta_train]
    )
    fwrite(
        x = dt,
        file = paste0(
            path2save,
            method_train[j, method],
            "_",
            method_train[j, meta_train],
            "_important_features.fasta"
        ),
        row.names = FALSE
    )
}
names(probe_lists) <- paste0(
    method_train[, method], "+", method_train[, meta_train]
)

# Convert probe IDs to gene entrez IDs using the feature_data file
feature_data <- fread("data/raw/GSE6477/feature_data.csv")
gene_lists_entrez <- lapply(probe_lists, function(li) {
    gene_ids <- unlist(
        strsplit(
            x = feature_data[rn %in% li, Gene.ID], split = "///",
            fixed = TRUE
        )
    )
    return(gene_ids)
})

# Perform enrichment analysis
go_results <- compareCluster(
    geneCluster = gene_lists_entrez,
    fun = "enrichGO",
    OrgDb = org.Hs.eg.db,
    ont = "BP", # Biological Process
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
)
fwrite(
    x = as.data.table(go_results),
    file = paste0(path2save, "enrichGO_results - GSE235356 vs allGLP96.csv")
)

# enrich pathways
pathway_results <- compareCluster(
    geneCluster = gene_lists_entrez,
    fun = ReactomePA::enrichPathway,
    organism = "human", # Homo sapiens
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
)
fwrite(
    x = as.data.table(pathway_results),
    file = paste0(
        path2save, "enrichPathway_results - GSE235356 vs allGLP96.csv"
    )
)

# KEGG pathway
kegg_results <- compareCluster(
    geneCluster = gene_lists_entrez,
    fun = "enrichKEGG",
    organism = "hsa", # Homo sapiens
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
)
fwrite(
    x = as.data.table(kegg_results),
    file = paste0(path2save, "enrichKEGG_results - GSE235356 vs allGLP96.csv")
)

# Disease Ontology Enrichment Analysis
do_results <- compareCluster(
    geneCluster = gene_lists_entrez,
    fun = "enrichDO",
    ont = "DO",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
)
fwrite(
    x = as.data.table(do_results),
    file = paste0(path2save, "enrichDO_results - GSE235356 vs allGLP96.csv")
)

# View results
terms_to_check <- paste0(
    "MAPK |RAS |RAF |MEK |ERK |PI3K |AKT |NF-KB |",
    "Jak-STAT|Wnt |Hedgehog|TNFa|mTOR",
    "|multiple myeloma",
    "|myeloid|leukemia|myeloma|Plasmacytoma|Amyloidosis",
    "|Chronic Lymphocytic Leukemia|POEMS|Heavy Chain Disease",
    "|Mastocytosis|Castleman|Lymphoma"
)

# Visualize GO Enrichment Analysis
dotplot(go_results) + ggtitle("GO Enrichment Analysis")
# Specific GO terms
go_results_df <- as.data.frame(go_results)
specific_go_results_df <- go_results_df[
    grepl(
        pattern = terms_to_check, x = go_results_df$Description,
        ignore.case = TRUE
    ),
]
specific_go_results <- new(
    "compareClusterResult",
    compareClusterResult = specific_go_results_df
)
go_plot <- dotplot(specific_go_results) +
    ggtitle("Selected GO Enrichment Analysis")
ggsave(
    filename = paste0(
        path2save,
        "Selected GO Enrichment Analysis - GSE235356 vs allGLP96.pdf"
    ),
    plot = go_plot, width = 15, height = 6
)

# Visualize Pathway Enrichment Analysis Results
pathway_results_df <- as.data.frame(pathway_results)
selected_pathway_results <- new(
    "compareClusterResult",
    compareClusterResult = pathway_results_df[
        grepl(
            pattern = terms_to_check,
            x = pathway_results@compareClusterResult$Description,
            ignore.case = TRUE
        ),
    ]
)
pathway_plot <- dotplot(selected_pathway_results) +
    ggtitle("Selected Reactome Pathways")
ggsave(
    filename = paste0(
        path2save,
        "Selected Pathway Enrichment Analysis - GSE235356 vs allGLP96.pdf"
    ),
    plot = pathway_plot, width = 12, height = 7
)

# Visualize KEGG Pathway Enrichment Results
kegg_df <- as.data.frame(kegg_results)
dotplot(kegg_results) + ggtitle("KEGG Pathway Enrichment Analysis")

specific_kegg_df <- kegg_df[
    grepl(
        pattern = terms_to_check,
        x = kegg_df$Description
    ),
]
specific_kegg_results <- new(
    "compareClusterResult",
    compareClusterResult = specific_kegg_df
)
kegg_plot <- dotplot(specific_kegg_results, showCategory = NULL) +
    ggtitle("Selected KEGG Pathway Enrichment Analysis")
ggsave(
    filename = paste0(
        path2save,
        "Selected KEGG Pathway Enrichment Analysis - GSE235356 vs allGLP96.pdf"
    ),
    plot = kegg_plot, width = 12, height = 4
)

# Visualize Disease Ontology Enrichment Results
dotplot(do_results) + ggtitle("Disease Ontology Enrichment Analysis")

# Focus on specific diseases
do_results_df <- as.data.frame(do_results)
inds <- grepl(
    pattern = terms_to_check,
    x = do_results_df$Description,
    ignore.case = TRUE
)
specific_do_results_df <- do_results_df[inds, ]
specific_do_results <- new(
    "compareClusterResult",
    compareClusterResult = specific_do_results_df
)
do_plot <- dotplot(specific_do_results) + ggtitle(
    "Selected Disease Ontology Enrichment Analysis"
)
ggsave(
    filename = paste0(
        path2save,
        "Selected Disease Ontology Enrichment Analysis - GSE235356 vs allGLP96.png"
    ),
    plot = do_plot, width = 16, height = 4
)
