import matplotlib.pyplot as plt
import os
from classes import Case
import numpy as np
from matplotlib import colors
import matplotlib.ticker as ticker


def colormap_tecplot_modern(num_levels=100):
    # tecplot modern colormap
    color_dict = {'red': [[0, 1.0, 1.0],
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

    colormap_modern = colors.LinearSegmentedColormap('ColormapTecplot', segmentdata=color_dict, N=num_levels)

    return colormap_modern


# plotting presets
def set_plot_globals():
    plt.rcParams['agg.path.chunksize'] = 10000
    plt.rcParams['font.size'] = 16


def plot_dc_bias_over_phase(xy, path_figures):
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


def plot_voltages_over_phase(xy_voltage1, xy_voltage2, xy_dc_bias, path_figures):
    path_save = os.path.join(path_figures, "Voltages")

    fig, ax = plt.subplots()
    ax.plot(xy_voltage1[0], xy_voltage1[1])
    ax.scatter(xy_voltage1[0], xy_voltage1[1])
    ax.plot(xy_voltage2[0], xy_voltage2[1])
    ax.scatter(xy_voltage2[0], xy_voltage2[1])
    ax.plot(xy_dc_bias[0], xy_dc_bias[1])
    ax.scatter(xy_dc_bias[0], xy_dc_bias[1])
    ax.set_xlim(0, 180)
    ax.set_xticks(xy_voltage1[0])
    ax.tick_params(direction='in', right=True, top=True)
    ax.set_xlabel("Phase angle (°)")
    ax.set_ylabel("Voltage (V)")
    ax.set_title("Voltages")
    plt.tight_layout()
    plt.savefig(path_save, dpi=600)
    plt.close()


def plot_ead(case, path_figures_iead, iead_max_energy=None, plot_species='ION-TOT'):
    """

    Args:
        case (Case): case object
        path_figures_iead (str): path to IEAD figure directory
        iead_max_energy (int): max energy of IEADs
        plot_species (str): name of species to plot

    Returns:

    """
    pos_plot_species = None
    for i, species in enumerate(case.pcmc_species):
        if species == plot_species:
            pos_plot_species = i
    if pos_plot_species is None:
        print(f"could not find species '{plot_species}'")
        return

    # plot ion tot
    array = case.pcmc_eads[pos_plot_species]

    plt.figure(figsize=(2, 5), dpi=600)
    ax = plt.axes()
    colormap = colormap_tecplot_modern()
    max_angle = case.pcmc_max_angle[0]
    max_energy = case.pcmc_max_energy[0]
    max_value = np.max(array)
    plt.imshow(array, aspect='auto', origin='lower', cmap=colormap,
               extent=[-max_angle, max_angle, 0, max_energy], norm=colors.LogNorm(max_value * 0.01, max_value))
    ax.xaxis.set_major_locator(ticker.LinearLocator(3))
    ax.xaxis.set_minor_locator(ticker.LinearLocator(5))
    ax.yaxis.set_major_locator(ticker.MaxNLocator(5))
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.tick_params(which="both", direction='in', right=True, top=True)
    ax.set_ylim(0, iead_max_energy)
    plt.tight_layout()
    file_name = "IEAD_" + plot_species + "_" + case.name
    path_file = os.path.join(path_figures_iead, file_name)
    plt.savefig(path_file)
