import pandas as pd

from fetch_data import fetch_train_test
from lasso.lasso_module import lasso_pipeline

# Input ####
filenames_expression = [
    "data/processed_gpl96_gpl570_affy44_platform/expression_rma.csv",
    "data/processed_gpl96_gpl570_affy44_platform/expression_binary_0.csv",
    "data/processed_gpl96_gpl570_affy44_platform/expression_binary_0.25.csv",
    "data/processed_gpl96_gpl570_affy44_platform/expression_binary_0.5.csv",
    "data/processed_gpl96_gpl570_affy44_platform/expression_binary_0.75.csv",
    "data/processed_gpl96_gpl570_affy44_platform/expression_ranking.csv",
    "data/processed_gpl96_gpl570_affy44_platform/expression_ratios.csv",
]
fmetadata_train = "data/processed_gpl96_gpl570_affy44_platform/metadata_train.csv"
fmetadata_holdout = "data/processed_gpl96_gpl570_affy44_platform/metadata_holdout.csv"

# Calculations ####
reps = 20
colnames_out = []
for cvi in range(5):
    colnames_out.append("cv." + str(cvi))
colnames_out.append("test")
colnames_out.append("filename")
file_expressioni = filenames_expression[0]
dfs_out = []
dfs_out_test_perdataset = []
for file_expressioni in filenames_expression:
    print(f"\n#### Reading file: {file_expressioni} ######")
    print("Fetch train test datasets")
    X_train, y_train, X_test, y_test = fetch_train_test(
        file_expression=file_expressioni,
        fmeta_train=fmetadata_train,
        fmeta_holdout=fmetadata_holdout,
    )
    print("Training Lasso")
    accuracies_cv_all = []
    report_cv_all = []
    accuracy_test_all = []
    report_test_all = []
    accuracy_test_per_dataset_all = []
    for repi in range(reps):
        print(f"- {repi} from {len(range(reps))}")
        (
            accuracies_cv,
            report_cv,
            accuracy_test,
            accuracy_test_per_dataset,
            report_test,
        ) = lasso_pipeline(
            X_train=X_train,
            X_test=X_test,
            y_train_in=y_train,
            y_test_in=y_test,
            kfold=5,
        )
        print(f"accuracies_cv: {accuracies_cv}")
        print(f"accuracy_test: {accuracy_test}")
        print(f"accuracy_test_per dataset: {accuracy_test_per_dataset}")
        print("report_test:")
        print(report_test)
        accuracies_cv_all.append(accuracies_cv)
        report_cv_all.append(report_cv)
        accuracy_test_all.append(accuracy_test)
        accuracy_test_per_dataset["rep"] = repi
        accuracy_test_per_dataset_all.append(accuracy_test_per_dataset)
        report_test_all.append(report_test)
    df_cv = pd.DataFrame(accuracies_cv_all)
    df_test = pd.DataFrame(accuracy_test_all)
    df_out = pd.concat([df_cv, df_test], axis=1)
    df_out["file"] = file_expressioni.split(sep="/")[2]
    df_out.columns = colnames_out
    dfs_out.append(df_out)
    df_out_test_pdt = pd.concat(accuracy_test_per_dataset_all, ignore_index=True)
    df_out_test_pdt["file"] = file_expressioni.split(sep="/")[2]
    dfs_out_test_perdataset.append(df_out_test_pdt)

df_all = pd.concat(dfs_out, ignore_index=True)
print("All performances")
print(df_all)
df_all.to_csv(
    "results/processed_gpl96_gpl570_affy44_platform/lasso_accuracy.csv", index=False
)

df_all_per_dt = pd.concat(dfs_out_test_perdataset, ignore_index=True)
print("All performances per dataset")
print(df_all_per_dt)
df_all_per_dt.to_csv(
    "results/processed_gpl96_gpl570_affy44_platform/lasso_accuracy_per_dataset.csv",
    index=False,
)
