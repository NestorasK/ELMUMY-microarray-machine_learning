import pandas as pd
import math
from collections import Counter
from itertools import combinations


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


def calculate_ratios(df_coli):
    # Generate all possible combinations of unique values
    value_combinations = list(combinations(df_coli, 2))
    # Compute the ratio for each combination
    my_ratios = []
    for a, b in value_combinations:
        my_ratios.append(a / b)
    return my_ratios


def calculate_ratios_df(df, verbose=True):
    """Calculate the ratios of all possible combinations per columns"""
    numberofcombinations = calc_combinations(n=df.shape[0], r=2)
    if verbose:
        print(f"It will calculate {numberofcombinations} ratios!")
    ratio_dfs = []
    for coli in df.columns:
        myratios = calculate_ratios(df_coli=df[coli])
        ratio_df = pd.DataFrame(myratios, columns=[coli])
        ratio_dfs.append(ratio_df)
    ratio_df_all = pd.concat(ratio_dfs, axis=1)
    if verbose:
        print("\nRatios dataframe:")
        print(ratio_df_all)
    return ratio_df_all


def calc_combinations(n, r):
    """Calculate number of combinations select 'r' from 'n'"""
    return math.factorial(n) // (math.factorial(r) * math.factorial(n - r))
