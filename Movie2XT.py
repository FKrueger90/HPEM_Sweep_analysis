import os, re
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy import ndimage

np.set_printoptions(threshold=100000)
np.set_printoptions(precision=3)


def movie2xt(name, path_movie1, path_save, freq, num_zones=0):
    """
    generates XT plot(s) of the Efield in z-direction.
    It uses finds tecplot movie files and generates the plots based on the "R"-averaged values

    Args:
        path_movie1 (str): path to movie.1.plt file
        path_save:
        num_zones:
        freq: fundamental frequency (MHz)
        name (str): case name
    """
    title = name + "_XT_Efield"

    # cells borders used for averaging area
    cell_boarder_bottom = 34
    cell_boarder_top = 50
    cell_boarder_left = 0
    cell_boarder_right = 80

    v_max = 50
    v_min = -v_max

    # find dimensions of Tecplot file
    # --------------------------------------------------------------------
    # find 'EZ', 'I' and 'J' variable list
    with open(path_movie1) as f:
        for line, row in enumerate(f):
            if 'EZ' in row:
                pos_in_zone = line - 2  # no of variable, R and Z omitted starting with 0
            # find I
            if re.search(r'.*I= *([0-9]*),', row):
                I = int(re.findall(r'.*I= *([0-9]*),', row)[0])
            # find J
            if re.search(r'.*J= *([0-9]*),', row):
                J = int(re.findall(r'.*J= *([0-9]*),', row)[0])
                break

    # Find max T if not provided
    if num_zones:
        Zones = num_zones
    else:
        # loop over lines beginning from last
        for line in reversed(list(open(path_movie1))):
            # fast coarse match
            if re.search(r'.*?T.*?', line):
                # precise match
                regexpr = r'.*=\s*(\d+).*$'
                if re.search(regexpr, line):
                    T = re.findall(regexpr, line)[0]
                    Zones = int(T)
                    break

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
                # print(f"ZONE {Zone}")
                line_zone = line + 2
                # print(f"line: {line_zone}" )

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
               extent=[0, 1 / freq, 0, 3.0], aspect=0.25)
    plt.colorbar()
    plt.xlabel(r"time ($\mu$s)")
    plt.ylabel('x (cm)')
    plt.tick_params(axis='both', which='both', labelcolor="black", tickdir='in', right=True)
    plt.title(title)
    y = np.linspace(0, 3.0, XT.shape[0])
    x = np.linspace(0, 1, XT.shape[1])
    contour = plt.contour(x, y, XT, colors='gray', levels=[0.0],
                          linestyles=['dashed', 'dashdot', 'dotted'], linewidths=1)
    plt.clabel(contour, inline=True, fontsize=8, fmt='%1.0f')

    if 'path_save' in locals():
        path_figure = path_save + f"\\{title}.png"
        if not os.path.isdir(path_save):
            os.mkdir(path_save)
        plt.savefig(path_figure, dpi=600)

        print(f"figure saved to {path_figure}\n")

    # plt.show()
    plt.close()
