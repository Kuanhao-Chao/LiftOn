import pycircos
import matplotlib
from matplotlib.colors import to_hex, to_rgb
import matplotlib.pyplot as plt
import numpy as np
import sys

target = sys.argv[1]
Garc = pycircos.Garc
Gcircle = pycircos.Gcircle

import collections

# Generate a color palette
chromosome_colors = plt.cm.tab20.colors  # Using matplotlib's tab20 colormap
chromosome_colors = chromosome_colors + plt.cm.Set3.colors
chromosome_color_dict = {}



# print("chromosome_colors: ", chromosome_colors)

N = 24
# test_cmaps = ['gist_rainbow','nipy_spectral','gist_ncar']
# segmented_cmaps = [matplotlib.colors.ListedColormap(plt.get_cmap(t)(np.linspace(0,1,N))) for t in test_cmaps]


# gradient = np.linspace(0, 1, 256)
# gradient = np.vstack((gradient, gradient))

# print("segmented_cmaps: ", segmented_cmaps)

# test = plt.get_cmap(segmented_cmaps[2])
# print(len(test))



# print(chromosome_colors.colors)

# # def plot_color_gradients(cmap_category, cmap_list, nrows):
# #     fig, axes = plt.subplots(nrows=nrows)
# #     fig.subplots_adjust(top=0.95, bottom=0.01, left=0.2, right=0.99)
# #     axes[0].set_title(cmap_category + ' colormaps', fontsize=14)

# #     for ax, name in zip(axes, cmap_list):
# #         ax.imshow(gradient, aspect='auto', cmap=plt.get_cmap(name))
# #         pos = list(ax.get_position().bounds)
# #         x_text = pos[0] - 0.01
# #         y_text = pos[1] + pos[3]/2.
# #         fig.text(x_text, y_text, name, va='center', ha='right', fontsize=10)

# #     # Turn off *all* ticks & spines, not just the ones with colormaps.
# #     for ax in axes:
# #         ax.set_axis_off()




# chromosome_colors = matplotlib.colors.ListedColormap(plt.get_cmap('gist_rainbow')(np.linspace(0,1,N))).colors

# chromosome_colors[:, :3] = 0.93 * chromosome_colors[:, :3]  # Adjust the saturation factor (0.7 in this case)
# chromosome_colors = np.clip(chromosome_colors, 0, 1)





# Assign colors to chromosomes
def assign_color_to_chromosome(chromosome):
    if chromosome not in chromosome_color_dict:
        color_index = len(chromosome_color_dict) % len(chromosome_colors)
        chromosome_color_dict[chromosome] = to_hex(chromosome_colors[color_index])

def blend_colors(color1, color2):
    """Blend two hex colors."""
    rgb1 = np.array(to_rgb(color1))
    rgb2 = np.array(to_rgb(color2))
    blended_rgb = (rgb1 + rgb2) / 2  # Average the RGB values
    return to_hex(blended_rgb)

# Function to add a curved arrow
def add_curved_arrow(circle, arc_id, start, end, raxis_range, color, direction):
    # Create a small arc to represent the curved arrow
    arrow_arc = Garc(arc_id=arc_id, size=end - start, interspace=0, raxis_range=raxis_range,
                     labelposition=1000, label_visible=False, facecolor=color)
    circle.add_garc(arrow_arc)

    # Add arrowhead as a label
    arrowhead_label = "→" if direction == "right" else "←"
    arrowhead_position = end if direction == "right" else start
    arrow_arc.add_label(text=arrowhead_label, position=arrowhead_position, 
                        labelposition=raxis_range[0] - 10, label_visible=True, labelsize=10)
    
# Set chromosomes and add labels

TYPES = ['ovp', 'nonovp']
for type in TYPES:
    circle = Gcircle(figsize=(8, 8))
    with open("/ccb/salz2/kh.chao/Lifton/experiments/miniprot_liftoff_circos_plot/human_refseq_chromosome_general_reverse_band.csv") as f:
        f.readline()
        for line in f:
            line = line.rstrip().split(",")
            chromosome_identifier = line[0]  # Full chromosome identifier (e.g., 1_r)
            name = chromosome_identifier.split('_')[0]  # Extract chromosome number
            assign_color_to_chromosome(name)  # Assign color
            length = int(line[-1])
            size_label = f"{round(length / 1e6)}M"  # Convert size to Mbp and format as a string
            label = f"{size_label}\n{name}"  # Concatenate size label and chromosome number (without _r or _t)

            # Add additional label for chromosome 1
            if chromosome_identifier == "1_r":
                label += "\nLiftoff ----->"
            elif chromosome_identifier == "1_t":
                label += "\n <----- miniprot"

            color = chromosome_color_dict[name]  # Get color for this chromosome

            # Main arc for the chromosome
            main_arc = Garc(arc_id=chromosome_identifier, size=length, interspace=1, raxis_range=(935, 985),
                            labelposition=70, label_visible=True, labelsize=7, facecolor=color, label=label)
            circle.add_garc(main_arc)

    # Set the positions of the arcs
    circle.set_garcs(0, 360)
    for arc_id in circle.garc_dict:
        circle.tickplot(arc_id, raxis_range=(985, 990), tickinterval=20000000, ticklabels=None)

    # Plot links
    # with open(f"/ccb/salz2/kh.chao/Lifton/results/{target}/visualization/circos/links.csv") as f:
    with open(f"/ccb/salz2/kh.chao/Lifton/results/{target}/visualization/circos/links_{type}.csv") as f:
        f.readline()
        for line in f:
            line = line.rstrip().split(",")
            name1 = line[0].split('_')[0]  # Extract chromosome number for the source
            name2 = line[3].split('_')[0]  # Extract chromosome number for the destination
            start1 = int(line[1]) - 1
            end1 = int(line[2])
            start2 = int(line[4]) - 1
            end2 = int(line[5])
            source = (line[0], start1, end1, 915)
            destination = (line[3], start2, end2, 915)

            if name1 not in chromosome_color_dict or name2 not in chromosome_color_dict:
                continue

            # Blend the colors of the source and destination chromosomes
            color1 = chromosome_color_dict[name1]
            color2 = chromosome_color_dict[name2]
            link_color = blend_colors(color1, color2)

            link_width = 0.5  # Adjust the line width as needed
            circle.chord_plot(source, destination, facecolor=link_color, edgecolor=link_color, linewidth=link_width)


    circle.figure

    # Save the figure to a file
    # plt.savefig(f'/ccb/salz2/kh.chao/Lifton/results/{target}/visualization/circos/circos.png', dpi=300, bbox_inches="tight")
    plt.savefig(f'/ccb/salz2/kh.chao/Lifton/results/{target}/visualization/circos/circos_{type}.png', dpi=300, bbox_inches="tight")
    plt.clf()
    plt.close()