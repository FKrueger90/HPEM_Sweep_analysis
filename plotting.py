import matplotlib.pyplot as plt
import os
from classes import Case
import numpy as np
from matplotlib import colors
import matplotlib.ticker as ticker

def colormap_tecplot_modern(num_levels=100):
    # tecplot modern colormap
    cdict = {'red': [[0, 1.0, 1.0],
                     [1 / 7, 1.0, 35.0 / 255],
                     [2 / 7, 0.0, 35.0 / 255],
                     [3 / 7, 0.0, 35.0 / 255],
                     [4 / 7, 0.0, 100.0 / 255],
                     [5 / 7, 1.0, 85.0 / 255],
                     [6 / 7, 217.0 / 255, 100.0 / 255],
                     [7 / 7, 1.0, 1.0]],
             'green': [[0, 35.0 / 255, 1.0],
                       [1 / 7, 0.0, 35.0 / 255],
                       [2 / 7, 0.0, 100.0 / 255],
                       [3 / 7, 1.0, 100.0 / 255],
                       [4 / 7, 1.0, 100.0 / 255],
                       [5 / 7, 1.0, 46.0 / 255],
                       [6 / 7, 117.0 / 255, 35.0 / 255],
                       [7 / 7, 0.0, 0.0]],
             'blue': [[0, 100.0 / 255, 1.0],
                      [1 / 7, 1.0, 100.0 / 255],
                      [2 / 7, 1.0, 100.0 / 255],
                      [3 / 7, 1.0, 35.0 / 255],
                      [4 / 7, 0.0, 35.0 / 255],
                      [5 / 7, 0.0, 10.0 / 255],
                      [6 / 7, 26.0 / 255, 35.0 / 255],
                      [7 / 7, 0.0, 0.0]]}

    cmap_modern = colors.LinearSegmentedColormap('testCmap', segmentdata=cdict, N=num_levels)

    return cmap_modern

# plotting presets
def set_plot_globals():
    plt.rcParams['agg.path.chunksize'] = 10000
    plt.rcParams['font.size'] = 16



def plot_dc_bias_over_phase(xy, path_figures, name=''):
    path_save = os.path.join(path_figures, "DCBias_Phase")
    x = xy[0]
    y = xy[1]
    fig, ax = plt.subplots()
    ax.plot(x, y)
    ax.scatter(x, y)
    ax.set_xlim(0, 180)
    ax.set_xticks(x)
    ax.tick_params(direction='in', right=True, top=True)
    ax.set_xlabel("Phase angle (°)")
    ax.set_ylabel("DC Self-Bias (V)")
    ax.set_title("DC Self-Bias")
    plt.tight_layout()
    plt.savefig(path_save, dpi=600)
    plt.close()


def plot_EAD(case):

    print("size:  ", len(case.pcmc_eads))
    # plot ion tot
    array = case.pcmc_eads[0]


    fig = plt.figure(figsize=(2, 5), dpi=600)
    ax = plt.axes()
    cmap = colormap_tecplot_modern()
    max_angle = case.pcmc_max_angle[0]
    max_energy = case.pcmc_max_energy[0]
    plt.imshow(array, aspect='auto', origin='lower', cmap=cmap, extent=[-max_angle, max_angle, 0, max_energy])
    ax.xaxis.set_major_locator(ticker.LinearLocator(3))
    ax.xaxis.set_minor_locator(ticker.LinearLocator(5))
    ax.yaxis.set_major_locator(ticker.MaxNLocator(5))
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.tick_params(which="both", direction='in', right=True, top=True)
    plt.tight_layout()
    plt.show()
