import pandas as pd
import matplotlib.pyplot as plt

# Read Expressions
df_expr = pd.read_csv("data/processed/geoSup_gpl96_platform.csv")
print("Expression values...")
print(df_expr)
df_expr.drop(["Unnamed: 0"], axis=1, inplace=True)
print(df_expr)

# Read binary file
thresholds = [0, 0.25, 0.5, 0.75]

for thresholdi in thresholds:
    filei = f"data/processed/geoSup_gpl96_platform_binary_{thresholdi}.csv"
    print(filei)
    df_binary = pd.read_csv(filei)
    print(df_binary)

    # Correlations
    corr_mat_expr = df_expr.drop(["rn"], axis=1).corr()
    corr_mat_binary = df_binary.drop(["ID"], axis=1).corr()

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
    plt.savefig(f"data/figures/cors_of_cors_{thresholdi}.pdf")
