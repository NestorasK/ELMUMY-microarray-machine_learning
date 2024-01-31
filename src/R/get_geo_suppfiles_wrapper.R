library(GEOquery)
library(affy)
library(oligo)
library(data.table)
get_geo_suppfiles_wrapper <- function(
    geo_dataset, folder2save,
    f_normexpression_cel) {
    cat("Download supplementary CEL files...\n")

    paths2supp_files <- getGEOSuppFiles(
        GEO = geo_dataset,
        makeDirectory = FALSE,
        baseDir = folder2save
    )

    dir.create(paste0(folder2save, "CEL_files"))
    untar(
        tarfile = rownames(paths2supp_files)[1],
        exdir = paste0(folder2save, "CEL_files")
    )

    cat("Calculate expression from CEL files\n")
    norm_expression <- tryCatch(
        {
            raw_data <- ReadAffy(
                celfile.path = paste0(folder2save, "CEL_files")
            )
            norm_expression <- affy::rma(
                raw_data,
                normalize = TRUE, background = TRUE
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
                normalize = TRUE, background = TRUE
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
        file = f_normexpression_cel
    )
}
