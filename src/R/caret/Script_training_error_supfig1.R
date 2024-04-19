rm(list = ls())
library(caret)
library(httpgd)
hgd()
hgd_browse()

path2read <- paste0(
    "results/experiments_caret/multiple_myeloma_stage/",
    "training_optimizing_multiclass_auc/"
)
models <- c("glmnet", "rf", "gbm", "svmLinear2", "svmRadial")

f_models <- vector(mode = "list", length = length(models))
counter <- 1
for (mli in models) {
    f_models[[counter]] <- list.files(
        path = path2read, pattern = paste0("^", mli, ".*\\.RData$"),
        full.names = TRUE
    )
    counter <- counter + 1
}
f_models <- unlist(f_models)
print(f_models)

models <- lapply(X = f_models, FUN = function(fi) {
    load(fi)
    return(modeli_fit)
})
names(models) <- sub(
    pattern = "_classes:['Normal', 'MGUS', 'MM']_dataset:['GSE6477']",
    replacement = "",
    x = sub(
        pattern = ".RData", replacement = "",
        x = basename(f_models)
    ),
    fixed = TRUE
)

resamps <- resamples(models)
models_comparison <- summary(resamps)
save(models_comparison,
    file = paste0(
        path2read, "all_models_comparison",
        ".RData"
    )
)
pdf(
    file = paste0(
        path2read, "supplementary_figure_1",
        ".pdf"
    ),
    width = 6,
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
        layout = c(1, 1),
        main = "Cross Validation performance"
    )
)
dev.off()
