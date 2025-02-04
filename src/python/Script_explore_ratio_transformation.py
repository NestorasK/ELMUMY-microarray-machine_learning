import pandas as pd
import matplotlib.pyplot as plt

# Read Expressions
df_expr = pd.read_csv("data/processed_gpl96_gpl570_affy44_platform/expression_rma.csv")
print("Expression values...")
print(df_expr)

# Read ratio file
df_ratio = pd.read_csv(
    "data/processed_gpl96_gpl570_affy44_platform/expression_ratiosfromranks.csv"
)
print("Ratio values...")
print(df_ratio)

# Correlations
print(
    "Calculate samples' correlation inside ratios from ranking and inside expression values"
)
print("Then calculate the correlations between the two correlation tables")
corr_mat_expr = df_expr.drop(["rn"], axis=1).corr()
corr_mat_ratio = df_ratio.corr()

corr_of_cors = []
for coli in corr_mat_expr.columns:
    corr_of_cors.append(corr_mat_expr[coli].corr(corr_mat_ratio[coli]))

plt.hist(corr_of_cors, color="lightgrey", edgecolor="black")
plt.title("Correlations between ratios from ranking and rma() expression values")
plt.xlabel("Pearson correlation")
plt.ylabel("Counts")
plt.savefig(
    f"data/processed_gpl96_gpl570_affy44_platform/cors_of_cors_ratiosfromranks.pdf"
)
