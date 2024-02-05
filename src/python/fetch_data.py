import pandas as pd


def fetch_train_test(file_expression, fmeta_train, fmeta_holdout):
    # Load expression
    df = pd.read_csv(file_expression)
    if "rn" in df.columns:
        df.index = df.rn
        df.drop("rn", axis=1, inplace=True)
    df = df.T
    df["rn"] = df.index
    # Make train
    meta_train = pd.read_csv(fmeta_train)
    df_train = pd.merge(left=df, right=meta_train[["rn", "class"]], on="rn")
    X_train = df_train.drop(["rn", "class"], axis=1)
    y_train = df_train["class"]
    # Make test
    meta_holdout = pd.read_csv(fmeta_holdout)
    df_test = pd.merge(left=df, right=meta_holdout, on="rn")
    X_test = df_test.drop(["rn", "class", "dataset"], axis=1)
    y_test = df_test[["class", "dataset"]]
    return (X_train, y_train, X_test, y_test)
