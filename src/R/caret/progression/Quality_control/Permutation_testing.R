rm(list = ls())
library(data.table)
library(ggplot2)
library(pROC)
library(parallel)
library(httpgd)
hgd()


path2save <- paste0(
    "results/experiments_caret/",
    "multiple_myeloma_progression/evaluation_optimizing_auc/"
)
path2read <- paste0(
    "results/experiments_caret/multiple_myeloma_progression/",
    "predictions_optimizing_auc/"
)
preds <- fread(file = paste0(path2read, "predictions.csv"))
predictions <- preds[
    dataset == "GSE235356" &
        meta_train != "metadata_train_classes:['MGUS', 'MM']_dataset:['GSE6477', 'EMTAB317'].csv",
]
predictions$training <- sub(
    pattern = "metadata_train_classes:['MGUS', 'MM']_dataset:[",
    replacement = "", x = predictions$meta_train, fixed = TRUE
)
predictions$training <- sub(
    pattern = "].csv",
    replacement = "", x = predictions$training, fixed = TRUE
)
predictions$training <- gsub(
    pattern = "'",
    replacement = "", x = predictions$training, fixed = TRUE
)
predictions$training <- gsub(
    pattern = ", ",
    replacement = " + ", x = predictions$training, fixed = TRUE
)


auc_holdout <- predictions[
    , .(
        auc_test = pROC::auc(
            response = class,
            predictor = MGUS
        )
    ),
    by = c("dataset", "method", "transformation", "training")
]

perm_process <- function(permi, predictions_copy) {
    predictions_copy$class <- predictions_copy$class[
        sample(x = 1:nrow(predictions_copy), size = nrow(predictions_copy))
    ]
    auc_holdout_permi <- predictions_copy[
        , .(
            auc_test = pROC::auc(
                response = class,
                predictor = MGUS
            )
        ),
        by = c("dataset", "method", "transformation", "training")
    ]
    colnames(auc_holdout_permi)[5] <- "auc_test_perm"
    return(auc_holdout_permi)
}
auc_holdout_perms <- mclapply(
    X = 1:1000, FUN = perm_process, predictions_copy = copy(predictions),
    mc.cores = 16
)
auc_perms <- do.call(
    what = cbind,
    args = lapply(auc_holdout_perms, FUN = function(dti) {
        return(dti[, 5])
    })
)
auc_holdout_truth_perms <- cbind(auc_holdout[, auc_test], auc_perms)
auc_holdout$pvalues <- apply(
    X = auc_holdout_truth_perms, MARGIN = 1,
    FUN = function(li) {
        sum(li[1] <= li) / length(li)
    }
)
set.seed(42)
plot_pvalues <- ggplot(
    data = auc_holdout,
    mapping = aes(x = transformation, y = pvalues, colour = training)
) +
    geom_jitter(height = 0, width = 0.2, shape = 1) +
    xlab("") +
    ylab("p-values") +
    facet_wrap(vars(method)) +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 60, hjust = 1),
        legend.position = "bottom"
    )


auc_holdout2plot <- cbind(auc_holdout[, -"pvalues"], auc_perms)
auc_holdout2plot <- melt.data.table(
    data = auc_holdout2plot,
    id.vars = c("dataset", "method", "transformation", "training")
)

plot_permdist <- ggplot(
    data = auc_holdout2plot[variable == "auc_test_perm", ],
    mapping = aes(
        x = transformation,
        y = value
    )
) +
    geom_jitter(
        height = 0, width = 0.2, shape = 1,
        colour = "grey"
    ) +
    geom_jitter(
        data = auc_holdout2plot[variable != "auc_test_perm", ],
        mapping = aes(
            x = transformation,
            y = value, colour = training
        ), width = 0.2, height = 0, shape = 1
    ) +
    ylab("AUC holdout") +
    xlab("") +
    facet_wrap(vars(method)) +
    theme_bw() +
    scale_y_continuous(breaks = seq(0.3, 0.9, 0.2), limits = c(0.3, 0.9)) +
    theme(
        axis.text.x = element_text(angle = 60, hjust = 1),
        legend.position = "bottom"
    )

plots <- cowplot::plot_grid(
    plotlist = list(plot_permdist, plot_pvalues), nrow = 2,
    labels = "AUTO"
)
ggsave(
    filename = paste0(path2save, "sanity_check_permutation_test.pdf"),
    plot = plots, width = 8, height = 11
)
