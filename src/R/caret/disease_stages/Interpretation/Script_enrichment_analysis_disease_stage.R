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
    "interpretation_optimizing_multiclass_auc/"
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
    myOverall > 0 & transformation != "ratios" & method != "svmLinear2",
]
importance_all <- importance_all[order(-myOverall), ]

probe_lists <- vector(
    mode = "list",
    length = length(unique(importance_all[, method]))
)
counter <- 1
for (methodi in unique(importance_all[, method])) {
    probe_lists[[counter]] <- unique(importance_all[method == methodi, rn])
    counter <- counter + 1

    # write file
    dt <- data.table(unique(importance_all[method == methodi, rn]))
    colnames(dt) <- paste0(">", methodi)
    fwrite(
        x = dt, file = paste0(path2save, methodi, "_important_features.fasta"),
        row.names = FALSE
    )
}
names(probe_lists) <- unique(importance_all[, method])

# Convert probe IDs to gene symbols using the annotation package
convert_probes_to_genes <- function(probes) {
    # Get the gene symbols
    gene_symbols <- unlist(mget(probes, hgu133plus2SYMBOL))
    return(gene_symbols)
}

# Apply the conversion function to all lists
gene_lists <- lapply(probe_lists, convert_probes_to_genes)

# Convert gene symbols to Entrez IDs
gene_lists_entrez <- lapply(gene_lists, function(genes) {
    bitr(genes,
        fromType = "SYMBOL",
        toType = "ENTREZID", OrgDb = org.Hs.eg.db
    )$ENTREZID
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

# KEGG pathway
kegg_results <- compareCluster(
    geneCluster = gene_lists_entrez,
    fun = "enrichKEGG",
    organism = "hsa", # Homo sapiens
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
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

# View results
View(as.data.frame(go_results))
View(as.data.frame(pathway_results))
View(as.data.frame(kegg_results))
View(as.data.frame(do_results))

# Visualize GO Enrichment Analysis
dotplot(go_results) + ggtitle("GO Enrichment Analysis")
# Specific GO terms
go_results_df <- as.data.frame(go_results)
specific_go_results_df <- go_results_df[
    grepl(pattern = "myeloid|leukemia", x = go_results_df$Description),
]
specific_go_results <- new(
    "compareClusterResult",
    compareClusterResult = specific_go_results_df
)
dotplot(specific_go_results) +
    ggtitle("Selected GO Enrichment Analysis")

# Visualize Pathway Enrichment Analysis Results
dotplot(pathway_results) + ggtitle("Pathway Enrichment Analysis")

# Visualize KEGG Pathway Enrichment Results
dotplot(kegg_results) + ggtitle("KEGG Pathway Enrichment Analysis")




# Visualize Disease Ontology Enrichment Results
dotplot(do_results) + ggtitle("Disease Ontology Enrichment Analysis")

# Focus on specific diseases
View(as.data.frame(do_results))
do_results_df <- as.data.frame(do_results)
inds <- grepl(
    pattern = "leukemia|myeloid|myeloma",
    x = do_results_df$Description
)
specific_do_results_df <- do_results_df[inds, ]
specific_do_results <- new(
    "compareClusterResult",
    compareClusterResult = specific_do_results_df
)
dotplot(specific_do_results, showCategory = NULL) + ggtitle(
    "Selected Disease Ontology Enrichment Analysis"
)
