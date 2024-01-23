import pandas as pd
import matplotlib.pyplot as plt

# Collect metadata

# "data/raw/GSE2113/phenodata.csv"
phenodata_GSE2113 = pd.read_csv("data/raw/GSE2113/phenodata.csv")
phenodata_GSE2113[["class", "other_info"]] = phenodata_GSE2113["title"].str.split(
    pat="-", regex=False, expand=True
)
phenodata_GSE2113["dataset"] = "GSE2113"
phenodata_GSE2113_clean = phenodata_GSE2113[["rn", "class", "dataset"]]


# "data/raw/GSE6477/phenodata.csv"
phenodata_GSE6477 = pd.read_csv("data/raw/GSE6477/phenodata.csv")
phenodata_GSE6477[["class_crude", "other_info"]] = phenodata_GSE6477["title"].str.split(
    pat="_", regex=False, expand=True
)
phenodata_GSE6477["dataset"] = "GSE6477"


def map_conditions(value):
    if value == "Relapsed MM":
        return "rMM"
    elif value == "MGUS":
        return "MGUS"
    elif value == "New MM":
        return "MM"
    elif value == "Smoldering MM":
        return "SMM"
    elif value == "rNew MM":
        return "MM"
    elif value == "Normal donor":
        return "Normal"
    else:
        return None


phenodata_GSE6477["class"] = phenodata_GSE6477["class_crude"].apply(map_conditions)
phenodata_GSE6477_clean = phenodata_GSE6477[["rn", "class", "dataset"]]
phenodata_GSE6477_clean


# "data/raw/GSE13591/phenodata.csv"
phenodata_GSE13591 = pd.read_csv("data/raw/GSE13591/phenodata.csv")
phenodata_GSE13591[["class_crude", "other_info"]] = phenodata_GSE13591[
    "title"
].str.split(pat="-", regex=False, expand=True)
phenodata_GSE13591["dataset"] = "GSE13591"
phenodata_GSE13591["class_crude"].unique()


def map_conditions(value):
    if value == "N":
        return "Normal"
    else:
        return value


phenodata_GSE13591["class"] = phenodata_GSE13591["class_crude"].apply(map_conditions)
phenodata_GSE13591_clean = phenodata_GSE13591[["rn", "class", "dataset"]]


# "data/raw/GSE14230/phenodata.csv"
phenodata_GSE14230 = pd.read_csv("data/raw/GSE14230/phenodata.csv")
phenodata_GSE14230[["class", "other_info"]] = phenodata_GSE14230["title"].str.split(
    pat="_", regex=False, expand=True
)
phenodata_GSE14230["dataset"] = "GSE14230"
phenodata_GSE14230_clean = phenodata_GSE14230[["rn", "class", "dataset"]]

phenodata_all = pd.concat(
    [
        phenodata_GSE13591_clean,
        phenodata_GSE14230_clean,
        phenodata_GSE2113_clean,
        phenodata_GSE6477_clean,
    ],
    ignore_index=True,
)
phenodata_all.to_csv("data/processed/metadata.csv", index=False)
class_freqs = phenodata_all["class"].value_counts()
class_freqs_perdataset = phenodata_all.groupby("dataset")["class"].value_counts()


# Plots
# Frequencies across all datasets
ax = class_freqs.plot(kind="bar", color="grey")
plt.title("Class counts")
plt.xlabel("Category")
plt.ylabel("Count")
plt.xticks(rotation=0, ha="center")
plt.savefig("data/figures/metadata_class_frequencies.pdf")

# Frequencies per dataset
grouped_counts_unstacked = class_freqs_perdataset.unstack()
# Plot the bar plot
ax = grouped_counts_unstacked.plot(kind="bar", stacked=True, colormap="viridis")

# Set labels and title
plt.xlabel("dataset")
plt.ylabel("Count")
plt.xticks(rotation=0, ha="center")
plt.title("Category per dataset")

# Save the plot
plt.savefig("data/figures/metadata_class_frequencies_perdataset.pdf")

# Add sums and save to csv
grouped_counts_unstacked["Sum"] = grouped_counts_unstacked.sum(axis=1)
grouped_counts_unstacked.loc["Sum"] = grouped_counts_unstacked.sum()
grouped_counts_unstacked

grouped_counts_unstacked.to_csv(
    "data/figures/metadata_class_frequencies_perdataset.csv"
)
