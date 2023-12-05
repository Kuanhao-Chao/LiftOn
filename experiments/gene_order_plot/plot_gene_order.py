from matplotlib import pyplot as plt
import numpy as np
import matplotlib as mpl
from matplotlib.axis import Axis 
import argparse
from liftofftools import filepaths
import warnings

COLOR_MAP = 'RdYlGn'
FIG_SIZE = [10, 10]
TICK_FONT_SIZE = 8
POINT_SIZE = 2
X_LABEL= 'Reference'
Y_LABEL = 'Target'
LABEL_SIZE = 12
GRID_WIDTH = 0.25
GRID_COLOR = 'gray'


def main(args):
    parser = argparse.ArgumentParser(description='plot gene synteny')
    parser.add_argument('input_file', help='tab seperated file with gene order output from synteny subcommand')
    args = parser.parse_args()
    order_arr = parse_input_file(args.input_file)
    print('Plotting gene order')
    plot_gene_order(order_arr, args)


def parse_input_file(input_file):
    order_arr = []
    with open(input_file) as f:
        for line in f:
            cols = line.strip().split("\t")
            order_arr.append([float(col) if col.replace(".", "", 1).isdigit() else col for col in cols])
    return order_arr


def plot_gene_order(order_arr, args, ref_chr_gene_ratio, tgt_chr_gene_ratio):
    print(">> plot_gene_order ...")
    output_file = filepaths.build_filepath([args.dir, filepaths.SYNTENY_OUTPUTS['plot']])
    if filepaths.make_file(output_file, args.force):
        x, y, c = get_scatter_points(order_arr)
        if len(x) == 0:
            mismatch_ids_warning = "No features with matching IDs to plot"
            warnings.warn(mismatch_ids_warning)
            return
        plt.rcParams["figure.figsize"] = FIG_SIZE
        mpl.rcParams['figure.dpi'] = 300
        color_boundaries = get_color_spacing()
        cmap_rb = plt.get_cmap(COLOR_MAP)
        colors = cmap_rb(np.linspace(0, 1, len(color_boundaries) - 1))
        cmap, norm = mpl.colors.from_levels_and_colors(color_boundaries, colors)
        fig, ax = plt.subplots()
        vlines, x_locs, x_labels = get_grid_and_ticks(order_arr, 1,  3)
        hlines, y_locs, y_labels = get_grid_and_ticks(order_arr, 2, 4)
        ax.vlines(x=vlines, ymin=0, ymax=np.max(y), color=GRID_COLOR, linestyle='--', linewidth=GRID_WIDTH)
        ax.hlines(y=hlines, xmin=0, xmax=np.max(x), color=GRID_COLOR, linestyle='--', linewidth=GRID_WIDTH)
        ax.set_axisbelow(True)
        cplot = ax.scatter(x, y, c=c, s=POINT_SIZE, cmap=cmap, norm=norm, zorder=5)
        cbar = fig.colorbar(cplot, ticks=color_boundaries, label='Protein Sequence Identity')
        cbar.ax.set_yticklabels([str(min(boundary, 1.0)) for boundary in color_boundaries])
        plt.xlim([-5, np.max(x) + 5])

        # Adjust x-axis labels to avoid overlapping
        plt.xlim([-5, np.max(x) + 5])
        plt.xticks(rotation=45, fontsize=TICK_FONT_SIZE, ha='right')


        # Set aspect ratio to be equal (make the figure square)
        ax.set_aspect('equal', adjustable='box')

        # Axis.set_major_locator(ax.xaxis, years)  

        # x_labels = [x_label for x_label in x_labels if len(x_label) > 5 else ""]
        x_labels = [x_label if ref_chr_gene_ratio[x_label] > 0.025 else "" for x_label in x_labels]
        y_labels = [y_label if tgt_chr_gene_ratio[y_label] > 0.025 else "" for y_label in y_labels]
        # ax.xaxis.set_major_locator(plt.FixedLocator(x_locs))
        # ax.xaxis.set_major_formatter(plt.FixedFormatter(x_labels))
        plt.xticks(x_locs, x_labels, rotation=45, fontsize=TICK_FONT_SIZE, ha='right')
        # ax.set_major_locator(x_locs, x_labels, rotation=45, fontsize=TICK_FONT_SIZE, ha='right')

        plt.ylim([-5, np.max(y) + 1])
        plt.yticks(y_locs, y_labels, fontsize=TICK_FONT_SIZE)
        plt.xlabel(X_LABEL, fontsize = LABEL_SIZE)
        plt.ylabel(Y_LABEL, fontsize = LABEL_SIZE)

        plt.tight_layout()
        plt.savefig(output_file, transparent=True)

        # plt.xticks(x_locs, x_labels, rotation=45, fontsize=TICK_FONT_SIZE, ha='right')

def get_scatter_points(order_arr):
    x,y,c = [],[],[]
    for row in order_arr:
        if row[1] != 'NA' and row[2] != 'NA':
            x.append(row[1])
            y.append(row[2])
            c.append(row[5])
    return x,y,c


def get_color_spacing():
    boundaries = 1 - np.linspace(0, 1, 10)
    # boundaries = 1 - np.logspace(-2.5, 0, 10)
    boundaries.sort()
    return np.around(boundaries,3)
    # return np.append(np.around(boundaries,3),1.001)


def get_grid_and_ticks(order_arr, sort_by_idx, seq_idx):
    subset_order_arr = [row for row in order_arr if row[sort_by_idx] != 'NA']
    subset_order_arr.sort(key= lambda x: x[sort_by_idx])
    previous_chrom = ""
    total_count = 0
    grid_lines,tick_locations,tick_labels = [], [], []
    for row in subset_order_arr:
        current_chrom = row[seq_idx]
        if current_chrom != previous_chrom:
            add_to_grid(grid_lines, total_count, tick_locations, tick_labels, previous_chrom)
        total_count += 1
        previous_chrom = current_chrom
    add_to_grid(grid_lines, total_count, tick_locations, tick_labels, previous_chrom)
    return grid_lines, tick_locations, tick_labels


def add_to_grid(grid_lines, total_count, tick_locs, tick_labels, chrom):
    grid_lines.append(total_count)
    if len(grid_lines) > 1:
        tick_locs.append(grid_lines[-1] - (grid_lines[-1] - grid_lines[-2]) / 2)
        tick_labels.append(chrom)


if __name__ == '__main__':
    main()