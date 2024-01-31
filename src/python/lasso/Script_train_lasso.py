import pandas as pd

from fetch_data import fetch_train_test
from lasso.lasso_module import lasso_pipeline

filenames_expression = [
    "data/processed/geoSup_gpl96_platform.csv",
    "data/processed/geoSup_gpl96_platform_binary_0.5.csv",
    "data/processed/geoSup_gpl96_platform_ranking.csv",
    "data/processed/geoSup_gpl96_platform_ratios.csv",
]

reps = 20
colnames_out = []
for cvi in range(5):
    colnames_out.append("cv." + str(cvi))
colnames_out.append("test")
colnames_out.append("filename")

dfs_out = []
for file_expressioni in filenames_expression:
    print(f"\n#### Reading file: {file_expressioni} ######")
    # print("Fetch train test")
    X_train, y_train, X_test, y_test = fetch_train_test(
        file_expression=file_expressioni,
        fmeta_train="data/processed/metadata_train.csv",
        fmeta_holdout="data/processed/metadata_holdout.csv",
    )
    print("Training Lasso")
    accuracies_cv_all = []
    accuracy_test_all = []
    accuracy_test_per_dataset_all = []
    for repi in range(reps):
        print(f"- {repi} from {len(range(reps))}")
        accuracies_cv, accuracy_test, accuracy_test_per_dataset = lasso_pipeline(
            X_train=X_train, X_test=X_test, y_train=y_train, y_test_in=y_test, kfold=5
        )
        print(f"accuracies_cv: {accuracies_cv}")
        print(f"accuracy_test: {accuracy_test}")
        print(f"accuracy_test_per dataset: {accuracy_test_per_dataset}")
        accuracies_cv_all.append(accuracies_cv)
        accuracy_test_all.append(accuracy_test)
        accuracy_test_per_dataset_all.append(accuracy_test_per_dataset)
    df_cv = pd.DataFrame(accuracies_cv_all)
    df_test = pd.DataFrame(accuracy_test_all)
    df_out = pd.concat([df_cv, df_test], axis=1)
    df_out["file"] = file_expressioni.split(sep="/")[2]
    df_out.columns = colnames_out
    dfs_out.append(df_out)

df_all = pd.concat(dfs_out, ignore_index=True)
print("All performances")
print(df_all)
df_all.to_csv("results/lasso/lasso_results.csv", index=False)
