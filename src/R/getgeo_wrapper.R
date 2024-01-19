library(GEOquery)
library(affy)
getgeo_wrapper <- function(geo_dataset, folder2save) {
    cat("Downloading geo provided expressions...\n")
    gset <- getGEO(
        GEO = geo_dataset,
        destdir = folder2save,
        GSEMatrix = TRUE, AnnotGPL = TRUE
    )

    if (length(gset) > 1) {
        idx <- grep("GPL6244", attr(gset, "names"))
    } else {
        idx <- 1
    }
    gset <- gset[[idx]]

    # make proper column names to match toptable
    fvarLabels(gset) <- make.names(fvarLabels(gset))

    cat("Extract and write...\n ")

    # Write feature data - probes to genes
    cat("- feature data...\n ")
    data.table::fwrite(
        x = data.table::data.table(gset@featureData@data,
            keep.rownames = TRUE
        ),
        file = paste0(folder2save, "feature_data.csv")
    )

    # Expression values
    cat("- expression data...\n ")
    # log2 transformation
    ex <- data.table::data.table(exprs(gset), keep.rownames = TRUE)
    colnames(ex)[1] <- "ID"
    data.table::fwrite(
        x = ex,
        file = paste0(folder2save, "expression_gset.csv")
    )

    # phenotypic data
    cat("- pheno data...\n ")
    phenodata <- data.table::data.table(
        pData(phenoData(gset)),
        keep.rownames = TRUE
    )
    data.table::fwrite(
        x = phenodata,
        file = paste0(folder2save, "phenodata.csv")
    )
}
