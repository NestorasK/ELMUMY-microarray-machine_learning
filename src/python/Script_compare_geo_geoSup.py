import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, pearsonr

# I selected the microarray that are under the same 'GPL96' platform
paths = [
    "data/raw/GSE2113/",
    "data/raw/GSE6477/",
    "data/raw/GSE13591/",
    "data/raw/GSE14230/",
]
path2savehists = "data/figures/"

for pathi in paths:
    print(f"\nworking on: {pathi}")
    f_gset_expression = pathi + "expression_gset.csv"
    f_gsetsup_expression = pathi + "norm_expression_from_cel.csv"

    gset_expression = pd.read_csv(filepath_or_buffer=f_gset_expression)
    gsetsup_expression = pd.read_csv(f_gsetsup_expression)

    # Histograms per sample - gset dataset
    # Plot a histogram with expression values
    gset_expression.hist(bins=10, figsize=(25, 25), edgecolor="black", grid=False)

    # Save the plot
    fi_hist = path2savehists + f_gset_expression.split(sep="/")[2] + "_geo_hist.pdf"
    print(f"saving histogram: '{fi_hist}'")
    plt.savefig(fi_hist, dpi=300)

    # Histograms per sample - gsetsup dataset
    # Plot a histogram with expression values
    gsetsup_expression.hist(bins=10, figsize=(25, 25), edgecolor="black", grid=False)

    # Save the plot
    fi_hist = path2savehists + f_gset_expression.split(sep="/")[2] + "_geosup_hist.pdf"
    print(f"saving histogram: '{fi_hist}'")
    plt.savefig(fi_hist, dpi=300)

    # Calculate per sample correlation between datasets
    pearsons = []
    pearsons_pvalues = []
    spearmans = []
    spearman_pvalues = []
    coli = 1
    for coli in range(1, gset_expression.shape[1]):
        coefi_p, pvaluei_p = pearsonr(
            x=gset_expression.iloc[:, coli], y=gsetsup_expression.iloc[:, coli]
        )
        pearsons.append(coefi_p)
        pearsons_pvalues.append(pvaluei_p)
        coefi_sp, pvaluei_sp = spearmanr(
            a=gset_expression.iloc[:, coli], b=gsetsup_expression.iloc[:, coli]
        )
        spearmans.append(coefi_sp)
        spearman_pvalues.append(pvaluei_sp)

    # plot Pearson
    plt.figure()
    plt.hist(x=pearsons, edgecolor="black", color="lightgrey")
    plt.xlabel("Pearson")
    plt.ylabel("Frequency")
    plt.title(
        "Samples Pearson correlation geo vs geosup expression\n" + pathi.split("/")[2]
    )
    filenamei = (
        path2savehists + pathi.split("/")[2] + "_pearson_geoVSgeosup_expression.pdf"
    )
    print(f"saving histogram pearson: {filenamei}")
    plt.savefig(filenamei)

    # plot Spearman
    plt.figure()
    plt.hist(x=spearmans, edgecolor="black", color="lightgrey")
    plt.xlabel("Spearman")
    plt.ylabel("Frequency")
    plt.title(
        "Samples Spearman correlation geo vs geosup expression\n" + pathi.split("/")[2]
    )
    filenamei = (
        path2savehists + pathi.split("/")[2] + "_spearman_geoVSgeosup_expression.pdf"
    )
    print(f"saving histogram sperman: {filenamei}")
    plt.savefig(filenamei)
