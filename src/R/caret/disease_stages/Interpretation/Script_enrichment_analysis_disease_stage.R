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
terms_to_check <- "MAPK|RAS|RAF|MEK|ERK|PI3K|AKT|NF-KB|STAT|Wnt|Hedgehog|TNFa|myeloid|leukemia|myeloma|cancer"

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
    ggtitle("Selected Pathway Enrichment Analysis")
ggsave(
    filename = paste0(path2save, "Selected Pathway Enrichment Analysis.pdf"),
    plot = pathway_plot, width = 7.5, height = 10
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
        path2save, "Selected KEGG Pathway Enrichment Analysis.pdf"
    ),
    plot = kegg_plot, width = 7.5, height = 9
)

# Visualize Disease Ontology Enrichment Results
dotplot(do_results) + ggtitle("Disease Ontology Enrichment Analysis")

# Focus on specific diseases
do_results_df <- as.data.frame(do_results)
inds <- grepl(
    pattern = terms_to_check,
    x = do_results_df$Description
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
    plot = do_plot, width = 7.5, height = 5
)
