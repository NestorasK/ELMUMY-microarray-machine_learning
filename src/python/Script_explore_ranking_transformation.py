import pandas as pd
import matplotlib.pyplot as plt

# Read Expressions
df_expr = pd.read_csv("data/processed_gpl96_gpl570_affy44_platform/expression_rma.csv")
print("Expression values...")
print(df_expr)

# Read ranking file
df_ranking = pd.read_csv(
    "data/processed_gpl96_gpl570_affy44_platform/expression_ranking.csv"
)
print("Ranking values...")
print(df_ranking)

# Correlations
print("Calculate samples' correlation inside ranking and inside expression values")
print("Then calculate the correlations between the two correlation tables")
corr_mat_expr = df_expr.drop(["rn"], axis=1).corr()
corr_mat_binary = df_ranking.drop(["rn"], axis=1).corr()

corr_of_cors = []
for coli in corr_mat_expr.columns:
    corr_of_cors.append(corr_mat_expr[coli].corr(corr_mat_binary[coli]))

plt.hist(corr_of_cors, color="lightgrey", edgecolor="black")
plt.title("Correlations between ranking and rma() expression values")
plt.xlabel("Pearson correlation")
plt.ylabel("Counts")
plt.savefig(f"data/processed_gpl96_gpl570_affy44_platform/cors_of_cors_ranking.pdf")
