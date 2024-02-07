from sklearn.model_selection import train_test_split, StratifiedKFold
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, classification_report
import pandas as pd


def lasso_pipeline(X_train, X_test, y_train_in, y_test_in, kfold):
    """Write something here"""
    # Create a Lasso model #####
    # Convert categorical labels to numerical encoding
    label_encoder = LabelEncoder()
    y_train = label_encoder.fit_transform(y_train_in)
    y_test = label_encoder.fit_transform(y_test_in["class"])
    # Split the data into training and testing sets
    # X_train, X_test, y_train, y_test = train_test_split(
    #     X, y_encoded, test_size=0.2, random_state=123, stratify=y_encoded
    # )
    # Standardize features (optional but recommended)
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)
    # Create a StratifiedKFold cross-validation strategy
    stratified_cv = StratifiedKFold(n_splits=kfold, shuffle=True)
    # Create and fit a Logistic Regression model with L1 regularization (Lasso)
    lasso_model = LogisticRegression(
        penalty="l1", solver="liblinear", multi_class="auto", max_iter=1000
    )
    # Perform stratified cross-validation
    accs_cv = []
    reports_cv = []
    for train_index, test_index in stratified_cv.split(X_train_scaled, y_train):
        X_train_cv, X_val_cv = X_train_scaled[train_index], X_train_scaled[test_index]
        y_train_cv, y_val_cv = y_train[train_index], y_train[test_index]
        lasso_model.fit(X_train_cv, y_train_cv)
        y_val_pred = lasso_model.predict(X_val_cv)
        accuracy = accuracy_score(y_val_cv, y_val_pred)
        accs_cv.append(accuracy)
        report = classification_report(y_val_cv, y_val_pred)
        reports_cv.append(report)
    # After cross-validation, you can train the final model on the entire training set
    lasso_model.fit(X_train_scaled, y_train)
    # Predict on the test set
    y_test_pred = lasso_model.predict(X_test_scaled)
    # Evaluate the performance on the test set
    accuracy_test = accuracy_score(y_test, y_test_pred)
    report_test = classification_report(y_test, y_test_pred)
    # accuracy per dataset
    accuracy_test_per_dataset = []
    for dti in y_test_in["dataset"].unique():
        inds = y_test_in["dataset"] == dti
        accuracy_test_per_dataset.append(
            accuracy_score(y_true=y_test[inds], y_pred=y_test_pred[inds])
        )
    accuracy_test_per_dataset_out = pd.DataFrame(
        {
            "dataset": list(y_test_in["dataset"].unique()),
            "accuracy_test": accuracy_test_per_dataset,
        },
    )
    return (
        accs_cv,
        reports_cv,
        accuracy_test,
        accuracy_test_per_dataset_out,
        report_test,
    )
