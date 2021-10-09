import os
import numpy as np
from Movie2XT import movie2xt
import plotting
from classes import Cases
#from classes import Case

np.set_printoptions(threshold=100000)

# config
# ------------------------------------------------------------------------
XT_colorbar = False  # show color-bars in XT-plots
XT_lower_half = False  # only show bottom half of XT plots

IEAD_max_energy = 3000  # max energy for IEAD in eV (None if not used)

# root directory
dir_root = os.path.abspath("D:\\UIGEL5_D_Florian\\Voltage_Waveform_Tailoring\\HPEM\\ArCF4O2"
                           "\\DarkSpaceGeometry\\2000_3000r1")

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

path_figures_XT_Field = os.path.join(path_figures, "XT_Field")
if not os.path.isdir(path_figures_XT_Field):
    os.mkdir(path_figures_XT_Field)

path_figures_IEAD = os.path.join(path_figures, "IEAD")
if not os.path.isdir(path_figures_IEAD):
    os.mkdir(path_figures_IEAD)

# set global plot parameters
plotting.set_plot_globals()

if not cases.constant_phase:
    # DC bias Plot
    xy_phase_bias = cases.get_value_pair("cwaveform_phase", "dc_bias", custom_waveform_only=True)
    plotting.plot_dc_bias_over_phase(xy_phase_bias, path_figures)
    # plot voltages
    xy_phase_voltage1 = cases.get_value_pair("cwaveform_phase", "final_voltages[0]", custom_waveform_only=True)
    xy_phase_voltage2 = cases.get_value_pair("cwaveform_phase", "final_voltages[1]", custom_waveform_only=True)
    plotting.plot_voltages_over_phase(xy_phase_voltage1, xy_phase_voltage2, xy_phase_bias, path_figures)

# load pcmc data
for case in cases:
    if case.path_pcmc is not None:
        case.read_pcmc_file()
        plotting.plot_ead(case, path_figures_IEAD, iead_max_energy=IEAD_max_energy, plot_species="ION-TOT")


# create XT-plots of E-field
print("\ncreating E-field XT plots...")
for case in cases:
    movie2xt(case, path_figures_XT_Field, lower_half=XT_lower_half, do_color_bar=XT_colorbar)
