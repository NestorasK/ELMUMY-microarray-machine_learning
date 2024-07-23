rm(list = ls())
library(GEOquery)
library(foreach)
library(data.table)
library(affy)

geo_datasets <- c(
    "GSE47552", "GSE2113", "GSE5900", "GSE2658", "GSE6477",
    "GSE80608", "GSE14230", "GSE24990", "GSE13591", "GSE235356",
    "GSE186537"
)

out <- foreach(geo_dataset = geo_datasets, .errorhandling = "pass") %do% {
    cat("\nWorking on", geo_dataset, "\n")

    folder2save <- paste0("data/raw/", geo_dataset, "/")
    f_normexpression_cel <- paste0(
        folder2save,
        "only_background_corrected_expression_from_cel.csv"
    )

    if (file.exists(f_normexpression_cel)) {
        cat("Skipping...\n")
    } else {
        dir.create(folder2save)
        cat("Folder to save:", folder2save, "\n")
        source("src/R/getgeo_wrapper.R")
        getgeo_wrapper(
            geo_dataset = geo_dataset,
            folder2save = folder2save
        )
        source("src/R/get_geo_suppfiles_wrapper.R")
        get_geo_suppfiles_wrapper(
            geo_dataset = geo_dataset,
            folder2save = folder2save,
            f_normexpression_cel = f_normexpression_cel
        )
    }
    return("all done")
}
