import os
import numpy as np
import plotting
from classes import Cases

np.set_printoptions(threshold=np.inf)
np.set_printoptions(linewidth=500)

# config
# ------------------------------------------------------------------------
# Energy angular distributions
plot_EADS = True  # plot particle EADS on surface
plot_EADS_species = ("ION-TOT", "E", "CF^")
IEAD_max_energy = 3000  # max energy for IEAD plot in eV (None if not used)
EEAD_max_energy = 100  # max energy for IEAD plot in eV (None if not used)

plot_EDFs = True  # integrated energy distribution functions
plot_ADFs = True  # integrated energy distribution functions

# electric field XT plots
plot_XT_Efield = True
XT_colorbar = False  # show color-bars in XT-plots
XT_lower_half = False  # only show bottom half of XT plots

# probe potential over time
plot_local_potential_over_time = True
potential_over_time_locations = [(40, 20), (40, 45), (40, 67)]
potential_over_time_labels = ["Bottom", "Bulk", "Top"]

# root directory
dir_root = os.path.abspath("D:\\UIGEL5_D_Florian\\Voltage_Waveform_Tailoring\\HPEM\\ArCF4O2"
                           "\\DarkSpaceGeometry\\2000_3000r2")

# script body
# ------------------------------------------------------------------------

# initialize case list object
cases = Cases()
# generate case objects by scanning the root folder
cases.scan_directory_for_cases(dir_root)
cases.determine_sweep()
cases.print_info()

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

path_figures_time_varying_potentials = os.path.join(path_figures, "time_varying_potentials")
if not os.path.isdir(path_figures_time_varying_potentials):
    os.mkdir(path_figures_time_varying_potentials)

# set global plot parameters
plotting.set_plot_globals()

if not cases.constant_phase:
    # plot DC bias over phase
    xy_phase_bias = cases.get_value_pair("cwaveform_phase", "dc_bias", custom_waveform_only=True)
    plotting.plot_dc_bias_over_phase(xy_phase_bias, path_figures)
    # plot voltage amplitudes over phase
    xy_phase_voltage1 = cases.get_value_pair("cwaveform_phase", "final_voltages[0]", custom_waveform_only=True)
    xy_phase_voltage2 = cases.get_value_pair("cwaveform_phase", "final_voltages[1]", custom_waveform_only=True)
    plotting.plot_voltages_over_phase(xy_phase_voltage1, xy_phase_voltage2, xy_phase_bias, path_figures)

# plot used geometry
plotting.plot_geometry(cases, path_figures,
                       plot_local_potential_over_time=plot_local_potential_over_time,
                       potential_over_time_locations=potential_over_time_locations)

# load particle energy-angular distributions from pcmc file and plot
if plot_EADS or plot_EDFs or plot_ADFs:
    cases.read_pcmc_file_all()

    # plot particle energy angular distributions
    if plot_EADS:
        for case in cases:
            for s in plot_EADS_species:
                plotting.plot_ead(case, path_figures_IEAD, iead_max_energy=IEAD_max_energy, plot_species=s)
    #if plot_EDFs:
    #    cases.gen
# plot comparative ion energy distributions
#for case in Cases:
#    case.geneate_EDFs()

# load and plot time varying voltages
if plot_local_potential_over_time:
    for case in cases:
        case.get_local_potential_over_time(potential_over_time_locations, potential_over_time_labels)
        plotting.plot_voltages_over_time(case, path_figures_time_varying_potentials)


# create XT-plots of E-field
if plot_XT_Efield:
    print("\ncreating E-field XT plots...")
    for case in cases:
        plotting.movie2xt(case, path_figures_XT_Field, lower_half=XT_lower_half, do_color_bar=XT_colorbar)
