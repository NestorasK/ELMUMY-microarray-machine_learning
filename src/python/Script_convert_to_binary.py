import pandas as pd
from transform_data import convert_to_binary_df, get_frequencies

# Convert expression values to binary
file_expression = "data/processed_glp96_gpl570_platform/expression_rma.csv"
print(f"\nWorking on file: {file_expression}")
df_expr = pd.read_csv(filepath_or_buffer=file_expression)
print("\nExpression values...")
print(df_expr)

print("\nConverting to binary...")
thresholds = [0, 0.25, 0.5, 0.75]

for thresholdi in thresholds:
    print(f"- threshold {thresholdi}")
    df_binary = convert_to_binary_df(
        df=df_expr.drop("rn", axis=1), quantile_threshold=thresholdi
    )
    df_binary.insert(loc=0, column="rn", value=df_expr["rn"])
    print(df_binary)

    freqs = get_frequencies(df=df_binary.drop(["rn"], axis=1))

    # Write binary file
    filenamei = (
        f"data/processed_glp96_gpl570_platform/expression_binary_{thresholdi}.csv"
    )
    print(f"Writing df_binary '{filenamei}'")
    df_binary.to_csv(path_or_buf=filenamei, index=False)
