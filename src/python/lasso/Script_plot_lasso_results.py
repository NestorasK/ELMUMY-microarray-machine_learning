import pandas as pd
import matplotlib.pyplot as plt

lasso_results = pd.read_csv("results/lasso_gpl96_platform/lasso_results.csv")
lasso_results["mean_cv_accuracy"] = lasso_results[
    ["cv.0", "cv.1", "cv.2", "cv.3", "cv.4"]
].mean(axis=1)
lasso_results["method"] = lasso_results["filename"].str.split(pat="_").str[3]
lasso_results["method"] = lasso_results["method"].str.split(pat=".").str[0]
lasso_results["method"][lasso_results["method"].isna()] = "normalized"
lasso_results["method"].unique()
lasso_results.rename(columns={"test": "test_accuracy"}, inplace=True)

fig, ax = plt.subplots()
lasso_results[["test_accuracy", "mean_cv_accuracy", "method"]].boxplot(
    by="method", grid=False, rot=45
)

plt.savefig("results/lasso_gpl96_platform/lasso_results.pdf")
