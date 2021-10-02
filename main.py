import os
import numpy as np
from Movie2XT import movie2xt
import plotting
from classes import Cases, Case

np.set_printoptions(threshold=100000)

# config
# ------------------------------------------------------------------------
XT_colorbar = False # show color-bars in XT-plots

# root directory
# dir_root = os.path.abspath("C:\\Users\\flori\\Desktop\\temp2\\ConstVoltage_200_1000")
#dir_root = os.path.abspath("D:\\UIGEL5_D_Florian\Voltage_Waveform_Tailoring\HPEM\ArCF4O2\DarkSpaceGeometry\ConstVoltage_200_1000")
dir_root = os.path.abspath("D:\\UIGEL5_D_Florian\\Voltage_Waveform_Tailoring\\HPEM\\ArCF4O2\\DarkSpaceGeometry\\1000_1000")
#dir_root = os.path.abspath("D:\\temp\\2000_2000")

# script body
# ------------------------------------------------------------------------

# initialize case list object
cases = Cases()
# generate case objects by scanning the root folder
cases.scan_directory_for_cases(dir_root)

print(f"{len(cases)} cases found:")
for case in cases:
    print("\n   " + case.name)
    print("      DC-bias:", case.dc_bias)
    print("      custom waveforms:", case.contains_custom)
    print("      phase:", case.cwaveform_phase)

cases.determine_sweep()

# Plotting
# ------------------------------------------------------------------------------
# generate figure folder if non existent
path_figures = os.path.join(dir_root, "Figures")
if not os.path.isdir(path_figures):
    os.mkdir(path_figures)

path_figures_XT_Field = os.path.join(path_figures, "_XT_Field")
if not os.path.isdir(path_figures_XT_Field):
    os.mkdir(path_figures_XT_Field)

# set global plot parameters
plotting.set_plot_globals()

if not cases.constant_phase:
    # DC bias Plot
    xy_phase_bias = cases.get_value_pair("cwaveform_phase", "dc_bias", custom_waveform_only=True)
    plotting.plot_dc_bias_over_phase(xy_phase_bias, path_figures)

# create XT-plots of Efield
print("\ncreating E-field XT plots...")
for case in cases:
    movie2xt(case, path_figures_XT_Field, lower_half=True, do_color_bar=XT_colorbar)
