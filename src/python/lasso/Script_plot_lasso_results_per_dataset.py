import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

lasso_results = pd.read_csv(
    "results/processed_gpl96_gpl570_affy44_platform/lasso_accuracy_per_dataset.csv"
)
lasso_results["transformation"] = (
    lasso_results["file"]
    .str.replace(pat="expression_", repl="", regex=False)
    .str.replace(pat=".csv", repl="", regex=False)
)
lasso_results["transformation"].value_counts()
lasso_results = lasso_results[lasso_results["transformation"] != "binary_0"]

# Create a FacetGrid with datasets as rows and separate plots for each dataset
dataset_order = lasso_results["dataset"].unique()
g = sns.FacetGrid(
    lasso_results,
    col="dataset",
    col_wrap=2,
    height=4,
    aspect=1.5,
    col_order=dataset_order,
)

# Map a scatter plot for each dataset
g.map(
    sns.stripplot,
    "transformation",
    "accuracy_test",
    order=lasso_results["transformation"].unique(),
    jitter=True,
    alpha=0.5,
)

# Set titles
g.set_titles("{col_name}")

# Set labels and title
g.set_axis_labels("Transformation", "Accuracy")
plt.subplots_adjust(top=0.9)

# Adjust layout to prevent overlapping labels
# plt.tight_layout()

# Set the suptitle outside the figure to prevent overlapping with the subplot titles
g.fig.suptitle("Accuracy vs Transformation for Each Dataset")

plt.savefig(
    "results/processed_gpl96_gpl570_affy44_platform/lasso_results_per_dataset_stripplot.pdf"
)
