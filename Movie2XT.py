import os
import numpy as np
import math
import matplotlib.pyplot as plt
import time
import sys
from scipy import ndimage
np.set_printoptions(threshold=100000)
np.set_printoptions(precision=3)
time_start = time.time()
# this file is used to generate one or a batch of XT plot(s) of the
# Efield in z-direction. It uses finds tecplot movie files and generates
# the plots based on the "R"-averaged values

# process parameters ( for plotting )
# -------------------------------------------------------------------------------

freq = 1  # (MHz)


# save to secondary loccation
path_secondary_save = r"D:\UIGEL5_D_Florian\Voltage_Waveform_Tailoring\Tecplot_movie1_to_XTplot\1000_1000out"


# finds all files containing 'keyword' in root and all subfolders and add to file_list
# -------------------------------------------------------------------------------
# path to root folder
root = r"D:\UIGEL5_D_Florian\Voltage_Waveform_Tailoring\HPEM\ArCF4O2\OldGeometry\1000W_top\1000_1000r1"
keyword = "movie1.pdt"
file_list = []
for path, subdirs, files in os.walk(root):
    for name in files:
        if keyword in name:
            file_list.append(os.path.join(path, name))
print(f"{len(file_list)} files found")


# main loop over all found files
# ==============================================================================
for file in file_list:
    # generate title
    print("file path:" + file)
    title = file[file.rfind('\\')+1:-11]
    print("title:" + title)

    # find 'EZ' in variable list
    with open(file) as f:
        for line, row in enumerate(f):
            if 'EZ' in row:
                PosInZone = line-2  # no of variable, R and Z omitted starting with 0

    I = 116
    J = 74
    Zones = 401
    LinesPerI = math.ceil(I/7)

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
    with open(file, "r") as fp:
        for line, rowtext in enumerate(fp):
            if 'ZONE' in rowtext:

                Zone+=1
                #print(f"ZONE {Zone}")
                line_zone = line+2
                #print(f"line: {line_zone}" )

                # because constants X and R are not repeated
                line_var_begin = line_zone + (J * LinesPerI * PosInZone) + (j * LinesPerI)
                if Zone == 0:
                    line_var_begin = line_zone + (J * LinesPerI * (PosInZone+2)) + (j * LinesPerI)
                line_var_end = line_var_begin + LinesPerI * J
                line_I_begin = line_var_begin
                line_I_end = line_I_begin + LinesPerI

            if line >= line_var_begin and line < line_var_end:

                if line>=line_I_begin and line<line_I_end:
                    row.extend([float(x) for x in rowtext.split()])

                if line == line_I_end-1:
                    # print(row)
                    array[j, :, Zone] = row
                    j+=1
                    line_I_begin = line
                    line_I_end = line_I_begin + LinesPerI+1
                    row = []
                if j == 74:
                    j = 0

    # create XT from 3D array
    cell_boarder_bottom = 34
    cell_boarder_top = 50
    cell_boarder_left = 0
    cell_boarder_right = 80
    print(array.shape)
    XT = np.mean(array[cell_boarder_bottom:cell_boarder_top, cell_boarder_left:cell_boarder_right, :], axis=1)
    splitpoint = 200
    XT = np.hstack((XT[:, splitpoint:], XT[:, :splitpoint]))
    # smooth in t-direction
    XT = ndimage.uniform_filter1d(XT, 20, 1)
    # plot heatmap
    plt.imshow(XT, cmap='jet', interpolation='bicubic', origin='lower', vmax=100, vmin=-100, extent=[0, 1/freq, 0, 3.0], aspect=0.25)
    plt.colorbar()
    plt.xlabel(r"time ($\mu$s)")
    plt.ylabel('x (cm)')
    plt.tick_params(axis='both', which='both', labelcolor="black", tickdir='in', right=True)
    plt.title(title)
    y = np.linspace(0, 3.0, XT.shape[0])
    x = np.linspace(0, 1, XT.shape[1])
    contour = plt.contour(x, y, XT, colors='gray', levels=[500.0, 5000, 10000.0],
                          linestyles=['dashed','dashdot','dotted'], linewidths=1)
    plt.clabel(contour, inline=True, fontsize=8, fmt='%1.0f')

    if 'path_secondary_save' in locals():
        path_figure = path_secondary_save + f"\\{title}.tiff"
        if not os.path.isdir(path_secondary_save):
            os.mkdir(path_secondary_save)
        plt.savefig(path_figure, dpi=600)

        print(f"figure saved to {path_figure}\n")


    #plt.show()
    plt.close()
