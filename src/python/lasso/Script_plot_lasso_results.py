import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

lasso_results = pd.read_csv("results/lasso_affy44_platform/lasso_accuracy.csv")
lasso_results["mean_cv_accuracy"] = lasso_results[
    ["cv.0", "cv.1", "cv.2", "cv.3", "cv.4"]
].mean(axis=1)
lasso_results["transformation"] = (
    lasso_results["filename"]
    .str.replace(pat="expression_", repl="", regex=False)
    .str.replace(pat=".csv", repl="", regex=False)
)
lasso_results["transformation"].value_counts()
lasso_results = lasso_results[lasso_results["transformation"] != "binary_0"]
lasso_results.rename(columns={"test": "test_accuracy"}, inplace=True)

lasso_melt = pd.melt(
    frame=lasso_results[["test_accuracy", "mean_cv_accuracy", "transformation"]],
    id_vars="transformation",
    var_name="metric",
    value_name="accuracy",
)

fig, ax = plt.subplots()
sns.boxplot(
    x="transformation",
    y="accuracy",
    hue="metric",
    hue_order=["mean_cv_accuracy", "test_accuracy"],
    data=lasso_melt,
    palette="Set2",
)

plt.title(label="Repeated 20 times stratified CV | lasso | MGUS, MM")
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig("results/lasso_affy44_platform/lasso_accuracy.pdf")
