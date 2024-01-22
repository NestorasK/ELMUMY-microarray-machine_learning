import pandas as pd
import matplotlib.pyplot as plt

# Read Expressions
df_expr = pd.read_csv("data/processed/geoSup_gpl96_platform.csv")
print("Expression values...")
print(df_expr)
df_expr.drop(["Unnamed: 0"], axis=1, inplace=True)
print(df_expr)

# Read ranking file
df_ranking = pd.read_csv("data/processed/geoSup_gpl96_platform_ranking.csv")
print(df_ranking)

# Correlations
corr_mat_expr = df_expr.drop(["rn"], axis=1).corr()
corr_mat_binary = df_ranking.drop(["ID"], axis=1).corr()

corr_of_cors = []
for coli in corr_mat_expr.columns:
    corr_of_cors.append(corr_mat_expr[coli].corr(corr_mat_binary[coli]))

plt.hist(corr_of_cors, color="lightgrey", edgecolor="black")
plt.title("Correlations between ranking and rma() expression values")
plt.xlabel("Pearson correlation")
plt.ylabel("Counts")
plt.savefig(f"data/figures/cors_of_cors_ranking.pdf")
