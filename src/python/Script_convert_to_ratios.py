import pandas as pd
import matplotlib.pyplot as plt
from transform_data import calculate_ratios_df
from file_operations import find_files


# Which genes to select
# - Keep the probs that have a stable ranking in normal samples
fmetadata = "data/processed_gpl96_gpl570_affy44_platform/metadata_train.csv"
print(f"metadata to read: {fmetadata}")
df_metadata = pd.read_csv(fmetadata)
df_metadata.rename({"rn": "samples"}, axis=1, inplace=True)
normal_samples = df_metadata[df_metadata["class"] == "Normal"]
print(f"Normal samples in training:")
print(normal_samples)

df_ranking = pd.read_csv(
    "data/processed_gpl96_gpl570_affy44_platform/"
    + "background_corrected_expression_ranking.csv"
)
df_ranking_normal_samples = df_ranking.loc[:, normal_samples["samples"]]
df_ranking_normal_samples["ranking_std"] = df_ranking_normal_samples.std(
    axis=1)
df_ranking_normal_samples.insert(loc=0, column="rn", value=df_ranking["rn"])
df_ranking_normal_samples.sort_values(by="ranking_std", inplace=True)

print("df_ranking_normal_samples:")
print(df_ranking_normal_samples)

# Select the most stable 210 probes.
# 210 probes so that the total features will be 21945 - close to the number
# of the total number of probes

# # Plots to test ###
# # Select 10 probes
# df_rank_select_normal = df_ranking_normal_samples.head(n=10)

# df_expr = pd.read_csv(
#     "data/processed_gpl96_gpl570_affy44_platform/"
#     + "background_corrected_expression_qnorm.csv"
# )

# df_expr_select_normal = df_expr.loc[
#     df_expr["rn"].isin(df_rank_select_normal["rn"]),
#     df_ranking_normal_samples.columns.drop("ranking_std"),
# ]
# transposed_df_expr_select_normal = df_expr_select_normal.drop(["rn"], axis=1).T
# transposed_df_expr_select_normal.boxplot()

# plt.title("Top 10 from the selected probes")
# plt.xlabel("Selected probes")
# plt.ylabel("Expression qnorm")
# plt.savefig(
#     "data/processed_gpl96_gpl570_affy44_platform/"
#     + "boxplots_selected_probes_ratios_trainingset.pdf"
# )

# # Unselected to test - 10 probes
# df_rank_unselect_normal = df_ranking_normal_samples.tail(n=10)
# df_expr_unselect_normal = df_expr.loc[
#     df_expr["rn"].isin(df_rank_unselect_normal["rn"]),
#     df_ranking_normal_samples.columns.drop("ranking_std"),
# ]
# transposed_df_expr_unselect_normal = df_expr_unselect_normal.drop(["rn"], axis=1).T

# plt.figure()
# transposed_df_expr_unselect_normal.boxplot()
# plt.title("Top 10 from the unselected probes")
# plt.xlabel("Probes")
# plt.ylabel("Expression qnorm")
# plt.savefig(
#     "data/processed_gpl96_gpl570_affy44_platform/"
#     + "boxplots_unselected_probes_ratios_trainingset.pdf"
# )

df_rank_select_normal = df_ranking_normal_samples.head(n=210)

# Convert expression values to ratios
file_expression = "data/processed_gpl96_gpl570_affy44_platform/" + \
    "background_corrected_expression_rma.csv"

print(f"\nReading file: {file_expression}")
df_expr = pd.read_csv(file_expression)
print("Expression values...")
print(df_expr)
# Convert expressions to ratios
df_expr_select = df_expr.loc[df_expr["rn"].isin(
    df_rank_select_normal["rn"]), :]
df_ratio = calculate_ratios_df(
    df=df_expr_select.drop(["rn"], axis=1), verbose=False
)
print(
    f"\nNumber of NAs values in ratio calculations: {
        df_ratio.isna().sum().sum()}"
)
print("\ndf_ratio:")
print(df_ratio)
filename_ratios = file_expression.replace(
    "background_corrected_expression_rma",
    "background_corrected_expression_ratios",
)
print(f"Writing file: {filename_ratios}")
df_ratio.to_csv(filename_ratios, index=False)
