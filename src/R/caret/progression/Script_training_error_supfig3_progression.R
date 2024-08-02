rm(list = ls())
library(caret)
library(dplyr)
library(httpgd)
hgd()
hgd_browse()

path2read_models <- paste0(
    "results/experiments_caret/multiple_myeloma_progression/",
    "training_optimizing_auc/"
)
models <- c("glmnet", "rf", "gbm", "svmLinear2", "svmRadial")

f_models <- vector(mode = "list", length = length(models))
counter <- 1
for (mli in models) {
    f_models[[counter]] <- list.files(
        path = path2read_models, pattern = paste0("^", mli, ".*\\.RData$"),
        full.names = TRUE
    )
    counter <- counter + 1
}
f_models <- unlist(f_models)
print(f_models)

models <- lapply(
    X = f_models,
    FUN = function(fi) {
        load(fi)
        return(modeli_fit)
    }
)
model_names <- sub(
    pattern = "classes:",
    replacement = "", x = basename(f_models),
    fixed = TRUE
)
# qnorm
model_names <- sub(
    pattern = "['MGUS', 'MM']_dataset:['EMTAB317']_['MGUS', 'MM']_dataset:['EMTAB317'].RData",
    replacement = "EMTAB317", x = model_names,
    fixed = TRUE
)
model_names <- sub(
    pattern = "_['MGUS', 'MM']_dataset:['GSE6477', 'EMTAB317']_['MGUS', 'MM']_dataset:['GSE6477', 'EMTAB317'].RData",
    replacement = "_GSE6477 + EMTAB317", x = model_names,
    fixed = TRUE
)
model_names <- sub(
    pattern = "['MGUS', 'MM']_dataset:['GSE6477']_['MGUS', 'MM']_dataset:['GSE6477'].RData",
    replacement = "GSE6477", x = model_names,
    fixed = TRUE
)
model_names <- sub(
    pattern = "_['MGUS', 'MM']_dataset:['GSE6477', 'GSE2113', 'EMTAB316', 'GSE13591']_['MGUS', 'MM']_dataset:['GSE6477', 'GSE2113', 'EMTAB316', 'GSE13591'].RData",
    replacement = "_GSE6477 + GSE2113 + EMTAB316 + GSE13591", x = model_names,
    fixed = TRUE
)
# other
model_names <- sub(
    pattern = "_['MGUS', 'MM']_dataset:['EMTAB317'].RData",
    replacement = "_EMTAB317", x = model_names,
    fixed = TRUE
)
model_names <- sub(
    pattern = "_['MGUS', 'MM']_dataset:['GSE6477', 'EMTAB317'].RData",
    replacement = "_GSE6477 + EMTAB317", x = model_names,
    fixed = TRUE
)
model_names <- sub(
    pattern = "_['MGUS', 'MM']_dataset:['GSE6477'].RData",
    replacement = "_GSE6477", x = model_names,
    fixed = TRUE
)
model_names <- sub(
    pattern = "_['MGUS', 'MM']_dataset:['GSE6477', 'GSE2113', 'EMTAB316', 'GSE13591'].RData",
    replacement = "_GSE6477 + GSE2113 + EMTAB316 + GSE13591", x = model_names,
    fixed = TRUE
)
names(models) <- model_names

# Calculate performance
train_dataset <- strsplit(
    x = model_names, split = "_", fixed = TRUE
) %>%
    sapply(
        FUN = function(li) {
            return(li[[length(li)]])
        }
    )

trdti <- "GSE6477 + GSE2113 + EMTAB316 + GSE13591"
for (trdti in unique(train_dataset)) {
    modelsi <- models[train_dataset == trdti]
    names(modelsi) <- sub(
        pattern = paste0("_", trdti),
        replacement = "", x = names(modelsi),
        fixed = TRUE
    )
    resamps <- resamples(modelsi)
    models_comparison <- summary(resamps)
    filenamei <- paste0(
        path2read_models, "all_models_comparison",
        "_with_train_dataset:", trdti,
        ".RData"
    )
    save(models_comparison, file = filenamei)
    pdf(
        file = sub(pattern = ".RData", ".pdf", filenamei),
        width = 9,
        height = 7
    )
    theme1 <- trellis.par.get()
    theme1$plot.symbol$col <- rgb(.2, .2, .2, .4)
    theme1$plot.symbol$pch <- 16
    theme1$plot.line$col <- rgb(1, 0, 0, .7)
    theme1$plot.line$lwd <- 2
    trellis.par.set(theme1)
    print(
        bwplot(
            resamps,
            layout = c(3, 1),
            main = paste("Cross Validation performance | training:", trdti)
        )
    )
    dev.off()
}
