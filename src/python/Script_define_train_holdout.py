# Based on the
# "data/figures/metadata_class_frequencies_perdataset.csv"
# table I decide to
# - use only the samples that belong in one of the following classes
#   - Normal
#   - MGUS
#   - MM
# - train using the samples in GSE6477 because the classes are more balanced
#   and test on the other datasets.

import pandas as pd

phenodata_all = pd.read_csv("data/processed_gpl96_gpl570_affy44_platform/metadata.csv")
data_all = phenodata_all[phenodata_all["class"].isin(["Normal", "MGUS", "MM"])]
train_all = data_all[data_all["dataset"] == "GSE6477"]
train_all["class"].value_counts()
hold_out = data_all[data_all["dataset"] != "GSE6477"]
hold_out[["dataset", "class"]].value_counts()
train_all.to_csv(
    "data/processed_gpl96_gpl570_affy44_platform/metadata_train.csv", index=False
)
hold_out.to_csv(
    "data/processed_gpl96_gpl570_affy44_platform/metadata_holdout.csv", index=False
)
