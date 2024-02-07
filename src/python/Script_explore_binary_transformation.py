import pandas as pd
import matplotlib.pyplot as plt

# Read Expressions
df_expr = pd.read_csv("data/processed_gpl96_gpl570_affy44_platform/expression_rma.csv")
print("Expression values...")
print(df_expr)

# Read binary file
thresholds = [0, 0.25, 0.5, 0.75]
for thresholdi in thresholds:
    filei = f"data/processed_gpl96_gpl570_affy44_platform/expression_binary_{thresholdi}.csv"
    print(filei)
    df_binary = pd.read_csv(filei)
    print(df_binary)

    # Correlations
    corr_mat_expr = df_expr.drop(["rn"], axis=1).corr()
    corr_mat_binary = df_binary.drop(["rn"], axis=1).corr()

    corr_of_cors = []
    for coli in corr_mat_expr.columns:
        corr_of_cors.append(corr_mat_expr[coli].corr(corr_mat_binary[coli]))

    plt.figure()
    plt.hist(corr_of_cors, color="lightgrey", edgecolor="black")
    plt.title(
        f"Correlations between binary and rma() expression values, thres:{thresholdi}"
    )
    plt.xlabel("Pearson correlation")
    plt.ylabel("Counts")
    plt.savefig(
        f"data/processed_gpl96_gpl570_affy44_platform/cors_of_cors_{thresholdi}.pdf"
    )
