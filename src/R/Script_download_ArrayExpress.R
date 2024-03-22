rm(list = ls())
library(ArrayExpress)
library(foreach)
library(data.table)
library(affy)
library(httpgd)
hgd()
hgd_browse()


# E-MTAB-316 ####
load("data/raw/E-MTAB-316/E-MTAB-316.eSet.RData")
# download manual from
# https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-316
expression_rma_emtab_316 <- rma(object = study[[1]], normalize = FALSE)

fwrite(
    x = data.table(
        expression_rma_emtab_316@assayData$exprs,
        keep.rownames = TRUE
    ),
    file = "data/raw/E-MTAB-316/only_background_corrected_expression_from_cel.csv"
)
fwrite(
    x = data.table(
        expression_rma_emtab_316@phenoData@data,
        keep.rownames = TRUE
    ),
    file = "data/raw/E-MTAB-316/phenodata.csv"
)


# E-MTAB-317 #####
# download manual from
# https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-317
rm(list = ls())
load("data/raw/E-MTAB-317/E-MTAB-317.eSet.RData")
expression_rma_emtab_317 <- rma(object = study, normalize = FALSE)

fwrite(
    x = data.table(
        expression_rma_emtab_317@assayData$exprs,
        keep.rownames = TRUE
    ),
    file = "data/raw/E-MTAB-317/only_background_corrected_expression_from_cel.csv"
)
fwrite(
    x = data.table(expression_rma_emtab_317@phenoData@data,
        keep.rownames = TRUE
    ),
    file = "data/raw/E-MTAB-317/phenodata.csv"
)
