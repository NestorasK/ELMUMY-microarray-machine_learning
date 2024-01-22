import pandas as pd
from transform_data import convert_to_binary_df, get_frequencies

# Convert expression values to binary
file_expression = "data/processed/geoSup_gpl96_platform.csv"
print(f"\nWorking on file: {file_expression}")
df_expr = pd.read_csv(filepath_or_buffer=file_expression)
print("Expression values...")
print(df_expr)
df_expr.drop(["Unnamed: 0"], axis=1, inplace=True)
print(df_expr)

print("Converting to binary...")
thresholds = [0, 0.25, 0.5, 0.75]

for thresholdi in thresholds:
    df_binary = convert_to_binary_df(
        df=df_expr.drop("rn", axis=1), quantile_threshold=thresholdi
    )
    df_binary.insert(loc=0, column="ID", value=df_expr["rn"])
    print(df_binary)

    freqs = get_frequencies(df=df_binary.drop(["ID"], axis=1))

    # Write binary file
    filenamei = f"data/processed/geoSup_gpl96_platform_binary_{thresholdi}.csv"
    print(f"Writing df_binary {filenamei}")
    df_binary.to_csv(path_or_buf=filenamei, index=False)
