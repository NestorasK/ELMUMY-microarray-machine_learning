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
classes_to_separate = ["MGUS", "MM"]
training_datasets = [["GSE6477"], ["EMTAB317"], ["GSE6477", "EMTAB317"]]
path2read = "data/processed_gpl96_gpl570_affy44_platform/"
path2write = "data/processed_gpl96_gpl570_affy44_platform/"

# Calculations
phenodata_all = pd.read_csv(path2read + "metadata.csv")
data_all = phenodata_all[phenodata_all["class"].isin(classes_to_separate)]
training_dataset = training_datasets[0]
for training_dataset in training_datasets:
    print(f"\ntraining_dataset:{training_dataset}")
    train_all = data_all[data_all["dataset"].isin(training_dataset)]
    train_all[["dataset", "class"]].value_counts()
    hold_out = data_all[~data_all["dataset"].isin(training_dataset)]
    hold_out[["dataset", "class"]].value_counts()
    filename_train = (
        path2write
        + "metadata_train_classes:"
        + str(classes_to_separate)
        + "_dataset:"
        + str(training_dataset)
        + ".csv"
    )
    print(f"filename_train: '{filename_train}'")
    train_all.to_csv(filename_train, index=False)
    filename_holdout = filename_train.replace("_train_", "_holdout_")
    print(f"filename_holdout: '{filename_holdout}'")
    hold_out.to_csv(filename_holdout, index=False)
