rm(list = ls())
library(GEOquery)
library(parallel)
library(data.table)
library(affy)

geo_datasets <- c(
    "GSE13591", "GSE14230", "GSE2113", "GSE6477",
    "GSE235356", "GSE5900"
    # "GSE47552",  "GSE2658",
    # "GSE80608",  "GSE24990",
    # "GSE186537"
)

myfun <- function(geo_dataset) {
    folder2save <- paste0("data/raw/", geo_dataset, "/")
    f_backgrexpression_cel <- paste0(
        folder2save, "background_corrected_expression_from_cel.csv"
    )

    cat("Calculate expression from CEL files\n")
    norm_expression <- tryCatch(
        {
            raw_data <- ReadAffy(
                celfile.path = paste0(folder2save, "CEL_files")
            )
            norm_expression <- affy::rma(
                raw_data,
                normalize = FALSE, background = TRUE
            )
        },
        error = function(err) {
            print(err)
            cat("Using package oligo\n")
            cel_files <- list.files(
                path = paste0(folder2save, "CEL_files/"), pattern = ".CEL.gz",
                full.names = TRUE
            )
            list.celfiles(
                path = paste0(folder2save, "CEL_files"), full.names = TRUE
            )
            raw_data <- read.celfiles(cel_files)
            norm_expression <- oligo::rma(
                raw_data,
                normalize = FALSE, background = TRUE
            )
            return(norm_expression)
        }
    )
    norm_expression_dt <- data.table(
        exprs(norm_expression),
        keep.rownames = TRUE
    )
    fwrite(
        x = norm_expression_dt,
        file = f_backgrexpression_cel
    )
}
mclapply(X = geo_datasets, FUN = myfun, mc.cores = 4)
