import pandas as pd
import matplotlib.pyplot as plt
from transform_data import calculate_ratios_df

# Convert expression values to ratios
file_expression = "data/processed/geoSup_gpl96_platform.csv"
print(f"Reading file: {file_expression}")
df_expr = pd.read_csv(file_expression)
print("Expression values...")
df_expr.drop(["Unnamed: 0"], axis=1, inplace=True)
df_expr.rename(columns={"rn": "ID"}, inplace=True)
print(df_expr)


# Which genes to select
# - Keep the probs that have a stable ranking in normal samples - TODO
df_metadata = pd.read_csv("data/processed/metadata.csv")
df_metadata["samples"] = df_metadata["rn"] + ".CEL.gz"
normal_samples = df_metadata[df_metadata["class"] == "Normal"]

df_ranking = pd.read_csv("data/processed/geoSup_gpl96_platform_ranking.csv")
df_ranking_normal_samples = df_ranking.loc[:, normal_samples["samples"]]
df_ranking_normal_samples["ranking_std"] = df_ranking_normal_samples.std(axis=1)
df_ranking_normal_samples.insert(loc=0, column="ID", value=df_ranking["ID"])
df_ranking_normal_samples.sort_values(by="ranking_std", inplace=True)

print("df_ranking_normal_samples:")
print(df_ranking_normal_samples)

# Select the most stable 210 probes.
# 210 probes so that the total features will be 21945 - close to the number
# of the total number of probes

# Plots to test ###
# Select 10 probes
df_rank_select_normal = df_ranking_normal_samples.head(n=10)
df_expr_select_normal = df_expr.loc[
    df_expr["ID"].isin(df_rank_select_normal["ID"]),
    df_ranking_normal_samples.columns.drop("ranking_std"),
]
transposed_df_expr_select_normal = df_expr_select_normal.drop(["ID"], axis=1).T
transposed_df_expr_select_normal.boxplot()

plt.title("Top 10 from the selected probes")
plt.xlabel("Selected probes")
plt.ylabel("Expression rma()")
plt.savefig("data/figures/boxplots_selected_probes_ratios.pdf")

# Unselected to test - 10 probes
df_rank_unselect_normal = df_ranking_normal_samples.tail(n=10)
df_expr_unselect_normal = df_expr.loc[
    df_expr["ID"].isin(df_rank_unselect_normal["ID"]),
    df_ranking_normal_samples.columns.drop("ranking_std"),
]
plt.figure()
transposed_df_expr_unselect_normal = df_expr_unselect_normal.drop(["ID"], axis=1).T
transposed_df_expr_unselect_normal.boxplot()

plt.title("Top 10 from the unselected probes")
plt.xlabel("Probes")
plt.ylabel("Expression rma()")
plt.savefig("data/figures/boxplots_unselected_probes_ratios.pdf")

# Convert expressions to ratios
df_rank_select_normal = df_ranking_normal_samples.head(n=210)
df_expr_select = df_expr.loc[df_expr["ID"].isin(df_rank_select_normal["ID"]), :]

df_ratio = calculate_ratios_df(df=df_expr_select.drop(["ID"], axis=1), verbose=False)
print("df_ratio:")
print(df_ratio)
df_ratio.to_csv("data/processed/geoSup_gpl96_platform_ratios.csv", index=False)
