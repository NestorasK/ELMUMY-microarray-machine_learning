import pandas as pd
import matplotlib.pyplot as plt
from functools import reduce

# I selected the microarray that are under the same 'GPL96' platform
filenames = [
    "data/raw/GSE2113/expression_data.csv",
    "data/raw/GSE6477/expression_data.csv",
    "data/raw/GSE13591/expression_data.csv",
    "data/raw/GSE14230/expression_data.csv",
]

path2savehists = "data/figures/"
dfs = []
for fi in filenames:
    dfi = pd.read_csv(filepath_or_buffer=fi)
    print(dfi.shape)
    dfs.append(dfi)


common_rows = reduce(lambda left, right: pd.merge(left, right, on="ID"), dfs)
print("Common Rows:")
print(common_rows)
print(common_rows.shape)
common_rows.to_csv(path_or_buf="data/processed/gpl96_platform.csv")
