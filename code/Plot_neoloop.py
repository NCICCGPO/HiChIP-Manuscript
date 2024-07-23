from neoloop.visualize.core import *
import cooler
plt.rcParams['svg.fonttype'] = 'none'
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from copy import copy
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

#-----------------
# Plot neoloops and neotads on reconstructed HiChIP matrix (Figure 5d)
#-----------------

palette = copy(plt.get_cmap('Reds'))
color_list = [mcolors.rgb2hex(palette(i)) for i in range(palette.N)]
last_color = color_list[len(color_list) -1]
color_list = [color_list[i] for i in range(0, 256,35)]
color_list.insert(9, last_color)
color_list.insert(0, "#ffffff")

with open('BRCA-8D1E6006-85CB-484A-8B5C-30766D90137B-X005-S02-B1-T1_H3K27ac_AA_assemblies.txt') as file:
    assemblies = [line.rstrip() for line in file]
    
clr = cooler.Cooler('BRCA-8D1E6006-85CB-484A-8B5C-30766D90137B-X005-S02-B1-T1_H3K27ac_10000.allValidPairs.cool')
assembly = [i for i in lines if 'Amplicon4_Cycle1_' in i][0]
vis = Triangle(clr, assembly, n_rows=6, figsize=(7,5), track_partition=[5, 0.3,0.4,0.6, 2, 0.5], correct='sweight')
vis.matrix_plot(vmin=0, colormap = LinearSegmentedColormap.from_list("",color_list), vmax = 0.01)
vis.plot_chromosome_bounds(linewidth=2)
vis.plot_loops('BRCA-8D1E6006-85CB-484A-8B5C-30766D90137B-X005-S02-B1-T1_H3K27ac.neo-loops_0.95.txt'
               , face_color='none', marker_size=40, cluster=True, filter_by_res=True)
vis.plot_neoTAD(color = "#1F78B4")
vis.plot_genes(filter_=['USP22','CDK12','ERBB2','PRKCA'],fontsize=9)
vis.plot_signal('H3K27ac 1D'
                , 'BRCA-8D1E6006-85CB-484A-8B5C-30766D90137B-X005-S02-B1-T1_H3K27ac_1D-signal.norm.bw'
                , label_size=10, data_range_size=8, color='#666666', max_value=45)

vis.plot_arcs(lw=1.5, cutoff='top', arc_color='#666666', gene_filter=["CDK12"]) 
vis.plot_chromosome_bar(name_size=11, coord_size=8, color_by_order=['#33A02C','#1F78B4', "#6A3D9A"])