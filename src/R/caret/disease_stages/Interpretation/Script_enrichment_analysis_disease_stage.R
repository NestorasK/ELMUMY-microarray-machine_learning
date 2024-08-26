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

path2read <- paste0(
    "results/experiments_caret/multiple_myeloma_stage/",
    "interpretation_optimizing_multiclass_auc/variable_importance/"
)

path2save <- paste0(
    "results/experiments_caret/multiple_myeloma_stage/",
    "interpretation_optimizing_multiclass_auc/enrichment_analysis/"
)

importance_all <- fread(
    file = paste0(
        path2read,
        "variable_importance_caret.csv"
    )
)
importance_all$myOverall <- rowSums(
    x = as.matrix(
        importance_all[, c("MGUS", "MGUS", "MM", "Normal", "Overall")]
    ),
    na.rm = TRUE
)
importance_all <- importance_all[
    myOverall > 0 & method != "svmLinear2",
]
importance_all <- importance_all[order(-myOverall), ]

probe_lists <- vector(
    mode = "list",
    length = length(unique(importance_all[, method]))
)
counter <- 1
for (methodi in unique(importance_all[, method])) {
    probe_lists[[counter]] <- unique(
        x = unlist(
            strsplit(
                x = importance_all[method == methodi, rn],
                split = "///", fixed = TRUE
            )
        )
    )
    counter <- counter + 1

    # write file
    dt <- data.table(
        unique(
            x = unlist(
                strsplit(
                    x = importance_all[method == methodi, rn],
                    split = "///", fixed = TRUE
                )
            )
        )
    )
    colnames(dt) <- paste0(">", methodi)
    fwrite(
        x = dt, file = paste0(path2save, methodi, "_important_features.fasta"),
        row.names = FALSE
    )
}
names(probe_lists) <- unique(importance_all[, method])

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
    file = paste0(path2save, "enrichGO_results.csv")
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
    file = paste0(path2save, "enrichPathway_results.csv")
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
    file = paste0(path2save, "enrichKEGG_results.csv")
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
    file = paste0(path2save, "enrichDO_results.csv")
)

# View results
terms_to_check <- paste0(
    "MAPK |RAS |RAF |MEK |ERK |ERK1 |ERK2 |PI3K |AKT |NF-KB |",
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
    filename = paste0(path2save, "Selected GO Enrichment Analysis.pdf"),
    plot = go_plot, width = 7.5, height = 5
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
    filename = paste0(path2save, "Selected Pathway Enrichment Analysis.pdf"),
    plot = pathway_plot, width = 7.5, height = 11
)

# Visualize KEGG Pathway Enrichment Results
kegg_df <- as.data.frame(kegg_results)
dotplot(kegg_results) + ggtitle("KEGG Pathway Enrichment Analysis")

specific_kegg_df <- kegg_df[
    grepl(
        pattern = terms_to_check,
        x = kegg_df$Description,
        ignore.case = TRUE
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
        path2save, "Selected KEGG Pathway Enrichment Analysis.pdf"
    ),
    plot = kegg_plot, width = 7.5, height = 5
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
        path2save, "Selected Disease Ontology Enrichment Analysis.pdf"
    ),
    plot = do_plot, width = 7.5, height = 4
)
