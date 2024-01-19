import pandas as pd
from collections import Counter


def convert_to_binary(value, threshold):
    """
    Convert the value to 0 or 1 if it is below or above threhold
    """
    if value > threshold:
        return 1
    else:
        return 0


def convert_to_binary_df(df, quantile_threshold):
    """
    Convert the values of a DataFrame to 0 and 1s.
    The 0s correspond to the minimum value per column
    """
    df_binary = df.apply(
        lambda col: col.apply(
            convert_to_binary, threshold=col.quantile(quantile_threshold)
        )
    )
    return df_binary


def get_frequencies(df, verbose=True):
    """Calculate the frequencies of each element"""
    element_frequencies = df.apply(lambda col: col.value_counts())
    if verbose:
        print("Element frequencies per column:")
        print(element_frequencies)
    return element_frequencies
