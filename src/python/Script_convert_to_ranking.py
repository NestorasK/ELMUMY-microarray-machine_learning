import pandas as pd

# Convert expression values to ranking
# Note: highest rank - highest expression
file_expression = (
    "data/processed_gpl96_gpl570_affy44_platform/"
    + "background_corrected_expression_rma.csv"
)
print(f"\nWorking on file: {file_expression}")
df_expr = pd.read_csv(filepath_or_buffer=file_expression)
print("Expression values...")
print(df_expr)

print("Ranking the expression values per sample")
df_ranking = df_expr.drop("rn", axis=1).rank()
df_ranking = df_ranking / df_ranking.shape[0]
df_ranking.insert(loc=0, column="rn", value=df_expr["rn"])
print(df_ranking)
print(df_ranking.describe())
df_ranking.to_csv(
    path_or_buf="data/processed_gpl96_gpl570_affy44_platform/"
    + "background_corrected_expression_ranking.csv",
    index=False,
)
