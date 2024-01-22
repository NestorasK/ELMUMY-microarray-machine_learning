import pandas as pd

# Convert expression values to ranking
file_expression = "data/processed/geoSup_gpl96_platform.csv"
print(f"\nWorking on file: {file_expression}")
df_expr = pd.read_csv(filepath_or_buffer=file_expression)
print("Expression values...")
print(df_expr)
df_expr.drop(["Unnamed: 0"], axis=1, inplace=True)
print(df_expr)

print("Ranking the expression values per sample")
df_ranking = df_expr.drop("rn", axis=1).rank()
df_ranking.insert(loc=0, column="ID", value=df_expr["rn"])
print(df_ranking)

df_ranking.to_csv("data/processed/geoSup_gpl96_platform_ranking.csv", index=False)
