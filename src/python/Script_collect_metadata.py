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


# "data/raw/GSE235356/phenodata.csv"
# Progression samples
phenodata_GSE235356 = pd.read_csv("data/raw/GSE235356/phenodata.csv")
phenodata_GSE235356["disease state:ch1"].value_counts()


def map_conditions(value):
    if value == "Stable MGUS":
        return "MGUS"
    elif value == "Progressing MGUS":
        return "progressing_MGUS"
    else:
        return value


phenodata_GSE235356["class"] = phenodata_GSE235356["disease state:ch1"].apply(
    map_conditions
)
phenodata_GSE235356["dataset"] = "GSE235356"
phenodata_GSE235356_clean = phenodata_GSE235356[["rn", "class", "dataset"]]

# "data/raw/GSE2658/phenodata.csv"
# Not sure if to use


# "data/raw/GSE5900/phenodata.csv"
phenodata_GSE5900 = pd.read_csv("data/raw/GSE5900/phenodata.csv")


def map_conditions_gse5900(value):
    if "MGUS" in value:
        return "MGUS"
    elif "smoldering" in value:
        return "SMM"
    elif "Healthy" in value:
        return "Normal"
    else:
        return "something is wrong"


phenodata_GSE5900["class"] = phenodata_GSE5900["title"].apply(map_conditions_gse5900)
phenodata_GSE5900["class"].value_counts()
phenodata_GSE5900["dataset"] = "GSE5900"
phenodata_GSE5900_clean = phenodata_GSE5900[["rn", "class", "dataset"]]


# data/raw/E-MTAB-316/phenodata.csv
phenodata_emtab316 = pd.read_csv("data/raw/E-MTAB-316/phenodata.csv")
phenodata_emtab316["Characteristics..DiseaseState."].unique()


def map_conditions_emtab316(value):
    if value == "multiple myeloma":
        return "MM"
    elif value == "monoclonal gammopathy of unknown significance (MGUS)":
        return "MGUS"
    else:
        return "Something is wrong"


phenodata_emtab316["class"] = phenodata_emtab316[
    "Characteristics..DiseaseState."
].apply(map_conditions_emtab316)
phenodata_emtab316["class"].value_counts()
phenodata_emtab316["dataset"] = "EMTAB316"
phenodata_emtab316_clean = phenodata_emtab316[["rn", "class", "dataset"]]


# data/raw/E-MTAB-317/phenodata.csv
phenodata_emtab317 = pd.read_csv("data/raw/E-MTAB-317/phenodata.csv")
phenodata_emtab317["Characteristics..DiseaseState."].value_counts()


def map_conditions_emtab317(value):
    if value == "myeloma":
        return "MM"
    elif value == "monoclonal gammopathy of unknown significance (MGUS)":
        return "MGUS"
    else:
        return "Something is wrong"


phenodata_emtab317["class"] = phenodata_emtab317[
    "Characteristics..DiseaseState."
].apply(map_conditions_emtab317)
phenodata_emtab317["class"].value_counts()
phenodata_emtab317["dataset"] = "EMTAB317"
phenodata_emtab317_clean = phenodata_emtab317[["rn", "class", "dataset"]]

# Combine metadata all
phenodata_all = pd.concat(
    [
        phenodata_GSE13591_clean,
        phenodata_GSE14230_clean,
        phenodata_GSE2113_clean,
        phenodata_GSE6477_clean,
        phenodata_GSE235356_clean,
        phenodata_GSE5900_clean,
        phenodata_emtab316_clean,
        phenodata_emtab317_clean,
    ],
    ignore_index=True,
)
phenodata_all.to_csv(
    "data/processed_gpl96_gpl570_affy44_platform/metadata.csv", index=False
)
class_freqs = phenodata_all["class"].value_counts()
class_freqs_perdataset = phenodata_all.groupby("dataset")["class"].value_counts()


# Plots
# Frequencies across all datasets
ax = class_freqs.plot(kind="bar", color="grey", rot=45)
plt.title("Class counts")
plt.xlabel("Category")
plt.ylabel("Count")
plt.tight_layout()
plt.savefig(
    "data/processed_gpl96_gpl570_affy44_platform/metadata_class_frequencies.pdf"
)

# Frequencies per dataset
grouped_counts_unstacked = class_freqs_perdataset.unstack()
# Plot the bar plot
ax = grouped_counts_unstacked.plot(
    kind="bar", stacked=True, colormap="viridis", figsize=(8, 6)
)

# Set labels and title
plt.xlabel("dataset")
plt.ylabel("Count")
plt.xticks(rotation=45, ha="center")
plt.title("Category per dataset")
plt.tight_layout()
# Save the plot
plt.savefig(
    "data/processed_gpl96_gpl570_affy44_platform/metadata_class_frequencies_perdataset.pdf"
)

# Add sums and save to csv
grouped_counts_unstacked["Sum"] = grouped_counts_unstacked.sum(axis=1)
grouped_counts_unstacked.loc["Sum"] = grouped_counts_unstacked.sum()
grouped_counts_unstacked

grouped_counts_unstacked.to_csv(
    "data/processed_gpl96_gpl570_affy44_platform/metadata_class_frequencies_perdataset.csv"
)
