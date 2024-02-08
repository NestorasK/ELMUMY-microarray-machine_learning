library(data.table)
fetch_data <- function(expres, meta_train, meta_holdout) {
    # Fetch data
    if ("rn" %in% colnames(expr)) {
        expr_t <- transpose(l = expr, keep.names = "rn", make.names = "rn")
    } else {
        expr_t <- transpose(l = expr, keep.names = "rn")
    }
    train_dt <- merge.data.table(x = meta_train, y = expr_t, by = "rn")
    holdout_dt <- merge.data.table(x = meta_holdout, y = expr_t, by = "rn")

    # make train x, y
    train_x <- as.matrix(train_dt[, -c("rn", "class", "dataset")])
    train_y <- train_dt[, class]

    # make hold out x, y
    holdout_x <- as.matrix(holdout_dt[, -c("rn", "class", "dataset")])
    holdout_y <- holdout_dt[, class]

    out <- list(
        train_x = train_x,
        train_y = train_y,
        holdout_x = holdout_x,
        holdout_y = holdout_y
    )
    return(out)
}
