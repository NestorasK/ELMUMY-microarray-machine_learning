# Load required packages
library(pROC) # For multiclass AUC calculation

# Define custom function to calculate multiclass AUC
multiclass_auc <- function(data, lev = NULL, model = NULL) {
    # Extract predicted probabilities for each class
    predictions <- as.matrix(
        data[, -c(1, 2)]
    ) # Exclude the first column (observed classes)

    # Extract observed classes
    observed_classes <- data$obs

    # Calculate AUC for each class against the rest
    class_auc <- sapply(levels(observed_classes), function(class) {
        # Create a binary outcome matrix for the current class
        binary_outcome <- ifelse(observed_classes == class, 1, 0)

        # Calculate ROC curve and AUC for the current class
        roc_obj <- pROC::roc(
            response = binary_outcome,
            predictor = predictions[, class], quiet = TRUE
        )
        auc(roc_obj)
    })

    # Calculate mean AUC across all classes
    mean_auc <- mean(class_auc)
    names(mean_auc) <- "multiclass_auc"
    return(mean_auc)
}
