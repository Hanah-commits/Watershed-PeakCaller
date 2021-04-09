import sys
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from Watershed_02 import *


def df_to_list(df, identifier):
    """
    :param df: dataframe with expression data of all genes
    :param identifier: name of gene
    :return: list containing expression data of interested gene
    """

    # extract expression data for each gene into a new dataframe
    # drop the column with the name of the gene
    list_gene = df.loc[df[0] == identifier].iloc[:, 1:].values.tolist()  # convert into a list of values
    list_gene = [val for sublist in list_gene for val in sublist]  # convert list of lists to a flat list

    return list(map(float, list_gene))


def plot_exp(gene_lists, legend_list):
    """
    :param gene_lists: list containing expression data of genes of interest
    :param legend_list: list of identifiers of gene
    :return: plot expression of genes
    """

    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    for gene in gene_lists:
        ax1.plot(gene)

    ax1.tick_params(labelbottom=False)
    ax1.xaxis.set_major_locator(MaxNLocator(10))


    ax1.set_title('Gene expression')
    ax1.set_xlabel('Pseudotime')
    ax1.set_ylabel('Expression')
    plt.gca().legend(legend_list)
    plt.savefig('Exp.png')


def plot_peak(gene_lists, legend_list):
    """
    :param gene_lists: list containing expression data of genes of interest
    :param legend_list: list of identifiers of gene
    :return: plot peaks
    """

    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    for gene in gene_lists:
        ax1.plot(gene)

    ax1.tick_params(labelbottom=False)
    ax1.xaxis.set_major_locator(MaxNLocator(10))


    ax1.set_title('Gene expression - Peaks')
    ax1.set_xlabel('Pseudotime')
    ax1.set_ylabel('Expression')
    plt.gca().legend(legend_list)
    plt.savefig('Peaks.png')





def main():
    """
    1. fill missing expression values in the data with 0
    2. extract expression data of each gene of interest into list
    3. peak-call using using Watershed on expression data of each gene
    4. Visualise before and after peak call
    """


    threshold = int(sys.argv[1])
    gene_ids = sys.argv[2:]

    # 1
    data = pd.read_csv('yeast_cell_cycle.txt', sep='\t', header=None)
    data = data.fillna(0)

    # 2
    gene_exp_data = []
    for gene_id in gene_ids:
        gene_exp_data.append(df_to_list(data, gene_id))

    # 3
    gene_exp_peaks = []
    for gene_exp in gene_exp_data:
        w1 = Watershed(gene_exp, threshold)
        w1.init_global()
        w1.watershed()
        gene_exp_peaks.append(w1.peak_calling())

    # 4
    plot_exp(gene_exp_data,  gene_ids)  # plot before peak-calling
    plot_peak(gene_exp_peaks, gene_ids)  # plot after peak-calling




if __name__ == "__main__":
    main()
