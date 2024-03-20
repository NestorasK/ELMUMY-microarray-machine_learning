rm(list = ls())
library(data.table)
library(glmnet)
library(pROC)
library(ggplot2)
library(httpgd)
hgd()
hgd_browse()

files_perf <- list.files(
    path = "results/processed_gpl96_gpl570_affy44_platform",
    pattern = "performance_all", full.names = TRUE
)

perfs <- vector(mode = "list", length = length(files_perf))
for (fi in seq_len(length(files_perf))) {
    perfi <- fread(files_perf[fi])
    perfi$training_dt <- unique(perfi[metric == "auc_cvmean", dataset])
    perfs[[fi]] <- perfi
}
perfs <- rbindlist(perfs)

plot_auc <- ggplot(
    data = perfs[metric %in% "AUC" & !dataset %in% c("EMTAB317", "GSE6477"), ],
    mapping = aes(x = transformations, y = performance, fill = training_dt)
) +
    geom_boxplot() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom"
    ) +
    facet_grid(cols = vars(dataset))
# scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1))

ggsave(
    filename = paste0(
        "results/processed_gpl96_gpl570_affy44_platform/",
        "boxplot_performance_AUC:['MGUS', 'MM'].pdf"
    ),
    plot = plot_auc, width = 18,
    height = 7
)
