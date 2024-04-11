rm(list = ls())
library(ggplot2)
library(data.table)
library(boot)
library(httpgd)
hgd()
hgd_browse()

predictions <- fread(
    file = paste0(
        "results/experiments_caret/multiple_myeloma_stage/",
        "predictions/predictions.csv"
    )
)

# Collect cv accuracy ####
path2read_models <- "results/experiments_caret/multiple_myeloma_stage/training/"
predictions$transformation[
    grepl(pattern = "qnorm", x = predictions$transformation)
] <-
    "qnorm_classes:['Normal', 'MGUS', 'MM']_dataset:['GSE6477']"
accuracy_cv <- vector(
    mode = "list",
    length = length(unique(predictions$transformation))
)
for (i in seq_len(length.out = length(unique(predictions$transformation)))) {
    transformationi <- unique(predictions$transformation)[i]
    load(
        file = paste0(
            path2read_models, "models_comparison_",
            transformationi, ".RData"
        )
    )
    acc_df <- data.table(models_comparison$statistics$Accuracy,
        keep.rownames = TRUE
    )
    acc_df$transformation <- transformationi
    accuracy_cv[[i]] <- acc_df
}
accuracy_cv_dt <- do.call(what = "rbind", args = accuracy_cv)
colnames(accuracy_cv_dt)[1] <- "method"

# Calculate test accuracy ####
accuracy_holdout <- predictions[
    , .(accuracy_test = mean(class == class_pred)),
    by = c("dataset", "method", "transformation")
]

# Compare ####
accuracy_all <- merge.data.table(
    x = accuracy_holdout, y = accuracy_cv_dt,
    by = c("method", "transformation")
)
accuracy_all[, diff_acccv_accholdout := Mean - accuracy_test]
accuracy_all$transformation[
    accuracy_all$transformation == "qnorm_classes:['Normal', 'MGUS', 'MM']_dataset:['GSE6477']"
] <- "qnorm"

ggplot(
    data = accuracy_all,
    mapping = aes(
        x = transformation, y = diff_acccv_accholdout,
        colour = method
    )
) +
    geom_jitter(width = 0.1, shape = 1) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom"
    ) +
    facet_grid(cols = vars(dataset))


ggplot(
    data = accuracy_all,
    mapping = aes(
        x = transformation, y = accuracy_test,
        colour = method
    )
) +
    geom_jitter(width = 0.1, shape = 1) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom"
    ) +
    facet_grid(cols = vars(dataset))
