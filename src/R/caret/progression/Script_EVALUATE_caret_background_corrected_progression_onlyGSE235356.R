rm(list = ls())
library(ggplot2)
library(data.table)
library(dplyr)
library(boot)
library(httpgd)
hgd()
hgd_browse()

# Read predictions ####
files_predictions <- list.files(
    path = paste0(
        "results/experiments_caret/multiple_myeloma_progression_onlyGSE235356/",
        "predictions_optimizing_auc/"
    ), full.names = TRUE
)
predictions <- rbindlist(
    l = lapply(X = files_predictions, FUN = fread)
)

# Path to save
path2save <- paste0(
    "results/experiments_caret/multiple_myeloma_progression_onlyGSE235356/",
    "evaluation_optimizing_auc/"
)

# Collect cv auc ####
path2read_models <- paste0(
    "results/experiments_caret/multiple_myeloma_progression_onlyGSE235356/",
    "training_optimizing_auc"
)
transf_ncv <- unique(predictions[, c("transformation", "ncv_fold")])

auc_cv <- vector(
    mode = "list",
    length = nrow(transf_ncv)
)

for (i in seq_len(length.out = length(auc_cv))) {
    transf_ncvi <- transf_ncv[i]
    transformationi <- transf_ncvi[, transformation]
    ncvi <- transf_ncvi[, ncv_fold]
    file_cv_perf <- paste0(
        path2read_models,
        "/models_comparison_",
        transformationi, "_rep", ncvi, ".RData"
    )
    load(file = file_cv_perf)
    auc_dt <- data.table(
        models_comparison$statistics$ROC,
        keep.rownames = TRUE
    )
    auc_dt$transformation <- transformationi
    auc_dt$ncv <- ncvi
    auc_cv[[i]] <- auc_dt
}
auc_cv_dt <- do.call(what = "rbind", args = auc_cv)
colnames(auc_cv_dt)[1] <- "method"
colnames(auc_cv_dt)[10] <- "ncv_fold"

# Calculate test auc ####
auc_holdout <- predictions[
    , .(
        auc_test = pROC::auc(
            response = class,
            predictor = MGUS
        )
    ),
    by = c("method", "transformation", "ncv_fold")
]

# Compare ####
auc_all <- merge.data.table(
    x = auc_holdout, y = auc_cv_dt,
    by = c("method", "transformation", "ncv_fold")
)
auc_all[
    ,
    diff_multiauccv_multiauc_holdout := Mean - auc_test
]
fwrite(x = auc_all, file = paste0(path2save, "auc_all.csv"))

# Performance ####
auc_all2plot <- melt.data.table(
    data = unique(
        x = auc_all[, c(
            "method", "transformation", "ncv_fold",
            "auc_test", "Mean"
        )]
    ), id.vars = c("method", "transformation", "ncv_fold"),
    value.name = "auc", variable.name = "metric"
)
auc_all2plot$metric <- as.character(auc_all2plot$metric)
auc_all2plot$metric[auc_all2plot$metric == "Mean"] <- "auc_cvmean"

plot_boxplot_perf <- ggplot(
    data = auc_all2plot,
    mapping = aes(
        x = transformation, y = auc,
        fill = metric
    )
) +
    geom_boxplot() +
    ylab("AUC") +
    theme(
        axis.text.x = element_text(angle = 60, hjust = 1),
        legend.position = "bottom"
    ) +
    facet_wrap(facets = vars(method), ncol = 2)

ggsave(
    filename = paste0(path2save, "boxplot_onlyGSE235356.png"),
    plot = plot_boxplot_perf, width = 6, height = 9
)
