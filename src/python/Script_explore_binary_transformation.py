import pandas as pd
from transform_data import convert_to_binary_df, get_frequencies
from file_operations import find_files

# files2read
files_expression = find_files(
    directory_path="data/raw/", filename_pattern="*expression_data.csv"
)

# Convert expression values to binary
for file_expression in files_expression:
    print(f"\nWorking on file: {file_expression}")
    df_expr = pd.read_csv(filepath_or_buffer=file_expression)

    print(f"nrows: {df_expr.shape[0]}")
    if df_expr.shape[0] == 0:
        print("Skipping...")
        continue

    print("Expression values...")
    print(df_expr)

    print("Converting to binary...")
    df_binary = convert_to_binary_df(
        df=df_expr.drop("ID", axis=1), quantile_threshold=0.5
    )
    df_binary.insert(loc=0, column="ID", value=df_expr["ID"])
    print(df_binary)

    freqs = get_frequencies(df=df_binary.drop(["ID"], axis=1))

    # Write binary file
    filenamei = "data/processed/binary/" + file_expression.split(sep="/")[2] + ".csv"
    print(f"Writing df_binary {filenamei}")
    df_binary.to_csv(path_or_buf=filenamei)
