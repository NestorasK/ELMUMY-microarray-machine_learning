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

# User input
classes_to_train = ["MGUS", "MM"]
classes_to_test = ["MGUS", "MM", "progressing_MGUS"]
training_datasets = [["GSE6477"], ["EMTAB317"], ["GSE6477", "EMTAB317"]]
path2read = "data/processed_gpl96_gpl570_affy44_platform/"
path2write = "data/processed_gpl96_gpl570_affy44_platform/"

# Calculations
phenodata_all = pd.read_csv(path2read + "metadata.csv")
for training_dataset in training_datasets:
    print(f"\ntraining_dataset:{training_dataset}")
    train_all = phenodata_all[phenodata_all["dataset"].isin(
        training_dataset) & phenodata_all["class"].isin(classes_to_train)]
    train_all[["dataset", "class"]].value_counts()
    hold_out = phenodata_all[~phenodata_all["dataset"].isin(training_dataset)]
    hold_out = phenodata_all[~phenodata_all["dataset"].isin(
        training_dataset) & phenodata_all["class"].isin(classes_to_test)]
    hold_out[["dataset", "class"]].value_counts()
    filename_train = (
        path2write
        + "metadata_train_classes:"
        + str(classes_to_train)
        + "_dataset:"
        + str(training_dataset)
        + ".csv"
    )
    print(f"filename_train: '{filename_train}'")
    train_all.to_csv(filename_train, index=False)
    filename_holdout = (
        path2write
        + "metadata_holdout_classes:"
        + str(classes_to_test)
        + "_dataset:"
        + str(training_dataset)
        + ".csv"
    )
    print(f"filename_holdout: '{filename_holdout}'")
    hold_out.to_csv(filename_holdout, index=False)
