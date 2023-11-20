
def main():
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import seaborn as sns
    import sys

    # Constants
    GRID_COLOR = 'grey'
    GRID_WIDTH = 0.5
    POINT_SIZE = 10
    LABEL_SIZE = 16
    TICK_FONT_SIZE = 14
    X_LABEL = 'Chromosome'
    Y_LABEL = 'Copy Number'

    # Input
    input_file = sys.argv[1]
    output_file = sys.argv[2]

    # Read in data
    df = pd.read_csv(input_file, sep='\t', header=None)
    df.columns = ['chrom', 'start', 'end', 'copy_number']
    df['chrom'] = df['chrom'].astype(str)
    df['copy_number'] = df['copy_number'].astype(float)
    df['chrom'] = df['chrom'].str.replace('chr', '')
    df['chrom'] = df['chrom'].replace('X', 23)
    df['chrom'] = df['chrom'].replace('Y', 24)
    df['chrom'] = df['chrom'].astype(int)
    df = df.sort_values(by=['chrom', 'start'])
    df['chrom'] = df['chrom'].astype(str)
    df['chrom'] = df['chrom'].str.replace('23', 'X')
    df['chrom'] = df['chrom'].str.replace('24', 'Y')

    # Plot
    fig, ax = plt.subplots()
    sns.scatterplot(x='chrom', y='copy_number', data=df, ax=ax, s=POINT_SIZE)
    plt.xlabel(X_LABEL, fontsize=LABEL_SIZE)
    plt.ylabel(Y_LABEL, fontsize=LABEL_SIZE)
    plt.xticks(fontsize=TICK_FONT_SIZE)
    plt.yticks(fontsize=TICK_FONT_SIZE)
    plt.savefig(output_file, transparent=True)

if __name__ == "__main__":
    main()