# Notes
# Seven from eleven microarray datasets were selected
# Platforms:
# - GPL96
# - GPL570
# The four were not selected
# - two errors in the download
# - two different probe ids

rm(list = ls())
library(data.table)

# GSE13591 ####
expression13591 <- fread(
    "data/raw/GSE13591/background_corrected_expression_from_cel.csv"
)
feature13591 <- fread("data/raw/GSE13591/feature_data.csv")
colnames(expression13591) <- sapply(
    X = strsplit(x = colnames(expression13591), split = ".", fixed = TRUE),
    FUN = function(li) li[1]
)

# GSE14230 ####
expression14230 <- fread(
    "data/raw/GSE14230/background_corrected_expression_from_cel.csv"
)
feature14230 <- fread("data/raw/GSE14230/feature_data.csv")
colnames(expression14230) <- sapply(
    X = strsplit(x = colnames(expression14230), split = ".", fixed = TRUE),
    FUN = function(li) li[1]
)

# GSE2113 ####
expression2113 <- fread(
    "data/raw/GSE2113/background_corrected_expression_from_cel.csv"
)
feature2113 <- fread("data/raw/GSE2113/feature_data.csv")
colnames(expression2113) <- sapply(
    X = strsplit(x = colnames(expression2113), split = ".", fixed = TRUE),
    FUN = function(li) li[1]
)

# GSE235356 ####
expression235356 <- fread(
    "data/raw/GSE235356/background_corrected_expression_from_cel.csv"
)
feature235356 <- fread("data/raw/GSE235356/feature_data.csv")
colnames(expression235356) <- sapply(
    X = strsplit(x = colnames(expression235356), split = "_", fixed = TRUE),
    FUN = function(li) li[1]
)

# GSE5900 #####
expression5900 <- fread(
    "data/raw/GSE5900/background_corrected_expression_from_cel.csv"
)
feature5900 <- fread("data/raw/GSE5900/feature_data.csv")
colnames(expression5900) <- sapply(
    X = strsplit(x = colnames(expression5900), split = "_", fixed = TRUE),
    FUN = function(li) li[1]
)


# GSE6477 ####
expression6477 <- fread(
    "data/raw/GSE6477/background_corrected_expression_from_cel.csv"
)
feature6477 <- fread("data/raw/GSE6477/feature_data.csv")
colnames(expression6477) <- sub(
    pattern = ".CEL.gz", replacement = "", x = colnames(expression6477)
)


# E-MTAB-316 ####
expression_emtab316 <- fread(
    "data/raw/E-MTAB-316/background_corrected_expression_from_cel.csv"
)


# E-MTAB-317 ####
expression_emtab317 <- fread(
    "data/raw/E-MTAB-317/background_corrected_expression_from_cel.csv"
)


# Merge expression datasets ####
merged_expression <- Reduce(
    function(...) merge(..., by = "rn", all = FALSE),
    list(
        expression13591, expression14230, expression2113,
        expression235356, expression5900, expression6477,
        expression_emtab316, expression_emtab317
    )
)

# Write ####
fwrite(
    x = merged_expression,
    file = paste0(
        "data/processed_gpl96_gpl570_affy44_platform/",
        "background_corrected_expression_rma.csv"
    )
)


# GSE186537 - not working ####
expression186537 <- fread("data/raw/GSE186537/norm_expression_from_cel.csv")
feature186537 <- fread("data/raw/GSE186537/feature_data.csv")
# expr_feat_gse186537 <- merge.data.table(
#     x = feature186537[, c("rn", "Gene.ID")], y = expression186537, by = "rn"
# )
# sum(expression186537$rn %in% expression13591$rn)

# GSE24990 not working - no ids in expression data! ####
expression24990 <- fread("data/raw/GSE24990/expression_gset.csv")
feature24990 <- fread("data/raw/GSE24990/feature_data.csv")

# GSE47552 - different probe ids #####
expression47552 <- fread("data/raw/GSE47552/norm_expression_from_cel.csv")
feature47552 <- fread("data/raw/GSE47552/feature_data.csv")
expr_feat_gse47552 <- merge.data.table(
    x = feature47552[, c("rn", "Gene.ID")], y = expression47552, by = "rn"
)

# GSE80608 - different probe IDs ####
expression80608 <- fread("data/raw/GSE80608/norm_expression_from_cel.csv")
feature80608 <- fread("data/raw/GSE80608/feature_data.csv")
expr_feat_gse80608 <- merge.data.table(
    x = feature80608[, c("rn", "GB_LIST")], y = expression80608, by = "rn"
)

# GSE2658 ####
# - Removed due to study design
expression2658 <- fread("data/raw/GSE2658/expression_gset.csv")
feature2658 <- fread("data/raw/GSE2658/feature_data.csv")
expression2658[, "ID"]
expr_feat_gse2658 <- merge.data.table(
    x = feature2658[, c("ID", "Gene.ID")],
    y = expression2658, by = "ID"
)
expr_feat_gse2658_log <- log2(expression2658[, -"ID"])
expr_feat_gse2658_log[, rn := expression2658$ID]
