from matplotlib import colors
import matplotlib.ticker as ticker

import os
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import ndimage
from classes import Case, Cases

np.set_printoptions(threshold=100000)
np.set_printoptions(precision=3)


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
    plt.style.use('classic')
    plt.rcParams['agg.path.chunksize'] = 10000
    plt.rcParams['font.size'] = 16
    plt.rcParams["xtick.direction"] = "in"
    plt.rcParams["ytick.direction"] = "in"
    plt.rcParams["xtick.top"] = True
    plt.rcParams["ytick.right"] = True
    # plot legend
    plt.rcParams['legend.fontsize'] = 'small'
    plt.rcParams['legend.loc'] = 'best'


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
    if plot_species == "E":
        file_name = "EEAD_" + case.name
    path_file = os.path.join(path_figures_iead, file_name)
    plt.savefig(path_file)
    plt.close()


def plot_voltages_over_time(case, path_figure):
    """
    generates voltage over time plots with electrode potential, bulk potential
    Args:
        case (Case): path to movie.1.plt file
        path_figure:
    """
    y_extension = 1.1  # extend y axis with buffer to avoid label collision at origin
    path_save = os.path.join(path_figure, "Potentials_over_Time" + case.name)
    fig, ax = plt.subplots()
    time_bins = case.potential_over_time[0].shape[0]
    time = np.linspace(0, 1e6/case.freq, time_bins )
    for i, potential in enumerate(case.potential_over_time):
        plt.plot(time, potential, label=case.potential_over_time_labels[i])

    y_min, y_max = ax.get_ylim()
    ax.set_ylim(y_min * y_extension, y_max * y_extension)
    ax.set_xlabel("Time (μs)")
    ax.set_ylabel("Potential (V)")
    ax.set_title("Time Varying Potentials")
    plt.tight_layout()
    plt.legend()
    plt.savefig(path_save, dpi=600)
    plt.close()


def plot_geometry(cases, path_figures, plot_local_potential_over_time=False, potential_over_time_locations=None):
    """
    checks if used mesh is identical in all cases. generates and saves all unique mesh figures
    Args:
        cases (Cases): Cases object
        path_figures (str): path to figure save directory
        plot_local_potential_over_time (bool): flag for adding potential probe points to plot
        potential_over_time_locations (list[tuple]): list of tuples defining the locations of probe points

    Returns:

    """
    # check if all meshes are identical
    meshes_identical = True
    for case in cases:
        if case.mesh is None:
            print("mesh not found\n aborting plotting")
            return None
        if not np.array_equal(case.mesh, cases[0].mesh):
            meshes_identical = False
    if meshes_identical:
        mesh = cases[0].mesh
        path_save = os.path.join(path_figures, "mesh")

        unique, rev = np.unique(mesh, return_inverse=True)
        rev = np.flip(rev.reshape(mesh.shape), axis=0)

        plt.figure()
        plt.imshow(rev, interpolation='none', cmap='gist_earth_r', origin='lower')
        if plot_local_potential_over_time:
            for loc in potential_over_time_locations:
                plt.scatter(loc[0], loc[1], s=200, marker='o', facecolors='w', edgecolors='k')
                plt.scatter(loc[0], loc[1], s=200, c='red', marker='x')

        plt.tight_layout()
        plt.ylim(0, rev.shape[0]-1)
        plt.xlim(0, rev.shape[1]-1)
        plt.savefig(path_save, dpi=600)
        plt.close()


def movie2xt(case, path_figure, lower_half=True, do_color_bar=True):
    """
    generates XT plot(s) of the Efield in z-direction.
    It uses finds tecplot movie files and generates the plots based on the "R"-averaged values

    Args:
        case (Case): path to movie.1.plt file
        path_figure:
    """
    title = case.name + "_XT_Efield"
    path_movie1 = case.path_movie1
    freq = case.freq*1e-6
    mesh = case.mesh
    cwafer = case.cwafer
    num_zones = case.rffac + 1

    gap = 3

    # find area for averaging
    # ---------------------------------------------------
    # default cell borders used for averaging area
    if mesh is None or cwafer is None:
        cell_boarder_bottom = 34
        cell_boarder_top = 50
        cell_boarder_left = 0
        cell_boarder_right = 80
    else:
        # find wafer width
        for num_line, line in enumerate(mesh[::-1]):
            if any(item in line for item in cwafer):
                cell_boarder_left = len(line)
                for num_cell, cell in enumerate(line):
                    if cell in cwafer:
                        if num_cell < cell_boarder_left:
                            cell_boarder_left = num_cell
                        cell_boarder_right = num_cell
                    else:
                        break

        # find electrode gap cells
        lines_bulk = []
        for num_line, line in enumerate(mesh[::-1]):
            if list(line[cell_boarder_left:cell_boarder_right]) == ['0'] * (cell_boarder_right - cell_boarder_left):
                lines_bulk.append(num_line)
        cell_boarder_top = max(lines_bulk)
        cell_boarder_bottom = min(lines_bulk)

        # only consider lower half
        if lower_half:
            cell_boarder_top = cell_boarder_bottom + int((cell_boarder_top - cell_boarder_bottom)/2)
            gap = gap/2

    v_max = 100
    v_min = -v_max

    # find dimensions of Tecplot file
    # --------------------------------------------------------------------
    # find 'EZ'
    with open(path_movie1) as f:
        for line, row in enumerate(f):
            if 'EZ' in row:
                pos_in_zone = line - 2  # no of variable, R and Z omitted starting with 0
                break

    I = case.movie_I
    J = case.movie_J
    Zones = case.movie_num_zones

    # read in data
    # --------------------------------------------------------------------
    LinesPerI = math.ceil(I / 7)

    # find lines of zones
    lines_zones = []
    row = []
    line_I = 0
    line_I_begin = 0
    line_I_end = 0
    array = np.zeros((J, I, Zones))
    Zone = -1
    j = 0
    line_var_begin = 0
    line_var_end = 0
    with open(path_movie1, "r") as fp:
        for line, rowtext in enumerate(fp):
            if 'ZONE' in rowtext:

                Zone += 1
                line_zone = line + 2

                # because constants X and R are not repeated
                line_var_begin = line_zone + (J * LinesPerI * pos_in_zone) + (j * LinesPerI)
                if Zone == 0:
                    line_var_begin = line_zone + (J * LinesPerI * (pos_in_zone + 2)) + (j * LinesPerI)
                line_var_end = line_var_begin + LinesPerI * J
                line_I_begin = line_var_begin
                line_I_end = line_I_begin + LinesPerI

            if line >= line_var_begin and line < line_var_end:

                if line >= line_I_begin and line < line_I_end:
                    row.extend([float(x) for x in rowtext.split()])

                if line == line_I_end - 1:
                    array[j, :, Zone] = row
                    j += 1
                    line_I_begin = line
                    line_I_end = line_I_begin + LinesPerI + 1
                    row = []
                if j == J:
                    j = 0

    # prepare data
    # --------------------------------------------------------------------
    # create XT from 3D array
    XT = np.mean(array[cell_boarder_bottom:cell_boarder_top, cell_boarder_left:cell_boarder_right, :], axis=1)
    splitpoint = 200
    XT = np.hstack((XT[:, splitpoint:], XT[:, :splitpoint]))
    # smooth in t-direction
    XT = ndimage.uniform_filter1d(XT, 20, 1)

    # plot heatmap
    # --------------------------------------------------------------------

    plt.imshow(XT, cmap='jet', interpolation='bicubic', origin='lower', vmax=v_max, vmin=v_min,
               extent=[0, 1 / freq, 0, gap], aspect="auto")
    plt.colorbar()
    plt.xlabel(r"time ($\mu$s)")
    plt.ylabel('x (cm)')
    plt.tick_params(axis='both', which='both', labelcolor="black", tickdir='in', right=True, top=True)
    plt.title(title)
    y = np.linspace(0, gap, XT.shape[0])
    x = np.linspace(0, 1, XT.shape[1])

    # add zero contour
    if (XT > 0).any():
        contour = plt.contour(x, y, XT, colors='gray', levels=[0.0],
                          linestyles=['dashed', 'dashdot', 'dotted'], linewidths=1)
        plt.clabel(contour, inline=True, fontsize=8, fmt='%1.0f')

    if 'path_figure' in locals():
        path_save = path_figure + f"\\{title}.png"
        if not os.path.isdir(path_figure):
            os.mkdir(path_figure)
        plt.savefig(path_save, dpi=600)
        print(f"   figure saved to {path_save}")
    plt.close()

    # plot with waveform:
    # ---------------------------------------------------------------------------------
    fig = plt.figure()
    # setup grid
    gs = gridspec.GridSpec(2, 2, height_ratios=[4, 1], width_ratios=[20, 1])
    gs.update(wspace=0.025, hspace=0.05)
    ax1 = plt.subplot(gs[0, 0])
    ax2 = plt.subplot(gs[1, 0], sharex=ax1)
    ax3 = plt.subplot(gs[:, 1])
    # ax setup
    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax2.get_yticklabels(), visible=False)
    ax1.tick_params(axis='both', which='both', labelcolor="black", tickdir='in', right=True, top=True)
    ax2.tick_params(axis='both', which='both', labelcolor="black", tickdir='in', right=True, top=True)
    ax2.tick_params(axis='y', right=False, left=False)

    # plot heatmap
    img = ax1.imshow(XT, cmap='jet', interpolation='bicubic', origin='lower', vmax=v_max, vmin=v_min,
                     extent=[0, 1 / freq, 0, gap], aspect='auto')

    # add zero contour
    if (XT > 0).any():
        contour = ax1.contour(x, y, XT, colors='gray', levels=[0.0],
                              linestyles=['dashed', 'dashdot', 'dotted'], linewidths=1)
        plt.clabel(contour, inline=True, fontsize=8, fmt='%1.0f')

    # generate waveform
    X = np.linspace(0, 1, 401)
    if case.contains_custom:
        Y = np.zeros(401)
        for i, harm in enumerate(case.custom_relharm):
            k = harm
            phi = case.custom_phase[i]
            Y += case.custom_relamp[i]*np.cos(2*k*np.pi*X + k*np.pi + phi/180*np.pi)
    else:
        Y = -np.sin(2*np.pi*X)
    ax2.plot(X, Y)

    # add color bar
    if do_color_bar:
        plt.colorbar(img, cax=ax3)
    else:
        ax3.axis("off")

    # axis labels
    ax2.set_xlabel(r"time ($\mu$s)")
    ax1.set_ylabel('x (cm)')
    plt.tick_params(axis='both', which='both', labelcolor="black", tickdir='in', right=True, top=True)
    fig.suptitle(title, fontsize=16)

    if 'path_figure' in locals():
        path_save = path_figure + f"\\{title}_waveform.png"
        if not os.path.isdir(path_figure):
            os.mkdir(path_figure)
        plt.savefig(path_save, dpi=600)
    plt.close()
    print(f"   figure saved to {path_save}")

