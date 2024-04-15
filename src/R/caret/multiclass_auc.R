# Load required packages
library(pROC) # For multiclass AUC calculation

# Define custom function to calculate multiclass AUC
multiclass_auc <- function(data, lev = NULL, model = NULL) {
    # Extract predicted probabilities for each class
    cols2select <- unique(as.character(data$obs))
    predictions <- as.matrix(
        data[, cols2select]
    )

    # Multiclass AUC
    multiauc <- auc(
        pROC::multiclass.roc(
            response = data$obs,
            predictor = predictions
        )
    )
    names(multiauc) <- "multiclass_auc"
    return(multiauc)
}
