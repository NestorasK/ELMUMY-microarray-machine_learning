import pandas as pd
import matplotlib.pyplot as plt
from transform_data import calculate_ratios_df

# Which genes to select
# - Keep the probs that have a stable ranking in normal samples
df_metadata = pd.read_csv("data/processed_gpl96_gpl570_affy44_platform/metadata.csv")
df_metadata.rename({"rn": "samples"}, axis=1, inplace=True)
normal_samples = df_metadata[df_metadata["class"] == "Normal"]

df_ranking = pd.read_csv(
    "data/processed_gpl96_gpl570_affy44_platform/expression_ranking.csv"
)
df_ranking_normal_samples = df_ranking.loc[:, normal_samples["samples"]]
df_ranking_normal_samples["ranking_std"] = df_ranking_normal_samples.std(axis=1)
df_ranking_normal_samples.insert(loc=0, column="rn", value=df_ranking["rn"])
df_ranking_normal_samples.sort_values(by="ranking_std", inplace=True)

print("df_ranking_normal_samples:")
print(df_ranking_normal_samples)

# Select the most stable 210 probes.
# 210 probes so that the total features will be 21945 - close to the number
# of the total number of probes

# Convert rankings to ratios
df_rank_select_normal = df_ranking_normal_samples.head(n=210)
df_ranking_select = df_ranking.loc[
    df_ranking["rn"].isin(df_rank_select_normal["rn"]), :
]
df_ratio = calculate_ratios_df(df=df_ranking_select.drop(["rn"], axis=1), verbose=False)

print(f"\nNumber of NAs values in ratio calculations: {df_ratio.isna().sum().sum()}")

print("\ndf_ratio:")
print(df_ratio)
df_ratio.to_csv(
    "data/processed_gpl96_gpl570_affy44_platform/expression_ratiosfromranks.csv",
    index=False,
)
