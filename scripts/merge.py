#############################################
########## 1. Import Libraries
#############################################
import sys
import pandas as pd
import numpy as np
import os
from scipy.sparse import csr_matrix
from scipy.io import mmwrite
import matplotlib.pyplot as plt
import seaborn as sns

#############################################
########## 2. Processing function
#############################################

## This function merges the different PCR batches, formats and saves the files
def merge_gene_count_files(count_files, output_basename):
    
    ### 1. Read counts and concatenate
    # Get files
    count_files_split = count_files.split(",")
    
    # Loop through files, count, and concatenate into a single dataframe
    print('Reading and concatenating count files...', flush=True)
    result_dataframe = pd.concat([pd.read_csv(file, sep='\t') for file in count_files_split])

    ### 2. Pivot and convert to sparse matrix
    # Gene gene_id or gene_name column
    gene_column = result_dataframe.columns[1]

    # Convert
    print('Converting to sparse matrix...', flush=True)
    pivot_dataframe = result_dataframe.pivot(index=gene_column, columns="barcode", values="unique_umi_count").fillna(0)

    # Convert to sparse matrix
    sparse_matrix = csr_matrix(pivot_dataframe.values)

    ### 3. Get cell and gene annotation
    # Get cell annotation
    cell_dataframe = pd.DataFrame({
        "barcode": pivot_dataframe.columns,
        "genes_detected": (sparse_matrix > 0).sum(axis=0).A1,  # Non-zero counts per cell,
        "total_counts": sparse_matrix.sum(axis=0).A1  # Total counts per cell
    })

    # Get cell annotation
    gene_dataframe = pd.DataFrame({
        gene_column: pivot_dataframe.index,
        "cells_expressed": (sparse_matrix > 0).sum(axis=1).A1  # Non-zero counts per cell
    })
    
    ### 4. Output
    print('Saving output files...', flush=True)
    # Get directories
    output_dir = os.path.dirname(output_basename)
    os.makedirs(output_dir, exist_ok=True)
    
    # Save as Matrix Market (MTX) format
    matrix_outfile = output_basename + "gene_counts.mtx"
    mmwrite(matrix_outfile, sparse_matrix)

    # Save metadata (cells)
    cell_annotation_file = output_basename + "cell_annotation.tsv"
    cell_dataframe.to_csv(cell_annotation_file, sep="\t", index=False)

    # Save metadata (genes)
    gene_annotation_file = output_basename + "gene_annotation.tsv"
    gene_dataframe.to_csv(gene_annotation_file, sep="\t", index=False)
    
    ### Create plot output for cells
    # Initialize plot dataframe
    plot_dataframe = cell_dataframe.sort_values('total_counts', ascending=False).reset_index(drop=True).reset_index()
    print(plot_dataframe.head())
    
    # Cells to highlight
    umi_highlight_threshold = 1000
    plot_dataframe['highlight'] = plot_dataframe['total_counts'] > umi_highlight_threshold
    
    # Define thresholds
    thresholds = [10000, 5000, 1000, 100]

    # Calculate the number of cells exceeding each threshold
    counts = [(plot_dataframe['total_counts'] > threshold).sum() for threshold in thresholds]

    # Generate dynamic annotation text
    text_annotation = "\n".join([f"Cells > {int(thresh):,} UMIs: {count:,}" for thresh, count in zip(thresholds, counts)])

    # Create the plot
    plt.figure(figsize=(7, 5))
    sns.lineplot(
        data=plot_dataframe,
        x="index",
        y="total_counts",
        hue="highlight",
        palette={True: "#1f78b4", False: "#a6cee3"},
        linewidth=1.5
    )

    # Set log scales
    plt.xscale('log')
    plt.yscale('log')

    # Customize y-axis ticks and add grid lines
    plt.grid(visible=True, which="major", axis="y", linestyle="--", linewidth=0.5, alpha=0.7)

    # Add labels, title, and legend
    plt.legend(title=f"UMIs > {int(umi_highlight_threshold):,}", loc="upper right")
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel("Barcodes", fontsize=12)
    plt.ylabel("UMI Count", fontsize=12)
    plt.title("Barcode Rank Plot", fontsize=14)

    # Add dynamic text annotation in the bottom-left corner
    plt.text(0.03, 0.05, text_annotation, fontsize=10, transform=plt.gca().transAxes, va='bottom', ha='left', bbox=dict(facecolor='white', edgecolor='grey', boxstyle='round,pad=0.5'))

    # Save to a PDF
    plt.tight_layout()
    barcode_plot_file = output_basename + "barcode_rank_plot.png"
    plt.savefig(barcode_plot_file, format="png", dpi=300)

    print('Output files saved successfully!', flush=True)

#############################################
########## 3. Main
#############################################

# Define main function
def main(args):
    
    # Run
    merge_gene_count_files(
        count_files = args.count_files,
        output_basename = args.output_basename
    )

#############################################
########## 4. Run
#############################################
# Run main function
if __name__ == "__main__":
    main(args)


# if __name__ == "__main__":
#     input_folder = sys.argv[1]
#     output_folder = sys.argv[2]
#     sampleID = sys.argv[3]
#     RT_matching_file = sys.argv[4]
#     merge_gene_count_files(input_folder, output_folder, sampleID, RT_matching_file)


