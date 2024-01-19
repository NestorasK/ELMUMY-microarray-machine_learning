import pandas as pd
from functools import reduce

# I selected the microarray that are under the same 'GPL96' platform
filenames = [
    "data/raw/GSE2113/norm_expression_from_cel.csv",
    "data/raw/GSE6477/norm_expression_from_cel.csv",
    "data/raw/GSE13591/norm_expression_from_cel.csv",
    "data/raw/GSE14230/norm_expression_from_cel.csv",
]

dfs = []
for fi in filenames:
    print(f"\nReading file: {fi}")
    dfi = pd.read_csv(filepath_or_buffer=fi)
    print(dfi)
    dfs.append(dfi)

common_rows = reduce(lambda left, right: pd.merge(left, right, on="rn"), dfs)
print("Common Rows:")
print(common_rows)
common_rows.to_csv(path_or_buf="data/processed/geoSup_gpl96_platform.csv")
