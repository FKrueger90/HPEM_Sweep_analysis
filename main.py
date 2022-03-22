# ╔═══════════════════════════════════════════════════════════════════════════════════╗ #
# ║   setup                                                                           ║ #
# ╚═══════════════════════════════════════════════════════════════════════════════════╝ #

# imports
import os
import numpy as np
import plotting
from classes import Cases
from classes import Config

# initial output setup
np.set_printoptions(threshold=np.inf)
np.set_printoptions(linewidth=500)

# ╔═══════════════════════════════════════════════════════════════════════════════════╗ #
# ║   config                                                                          ║ #
# ╚═══════════════════════════════════════════════════════════════════════════════════╝ #

# Energy angular distributions
plot_EADS = True                             # plot particle EADs on surface
plot_EDFs = True                             # plot integrated energy distribution functions
plot_ADFs = True                             # plot integrated energy distribution functions
plot_EADS_species = ("ION-TOT", "E", "CF^")  # pcmc species to be plotted
IEAD_max_energy = 3000                       # max energy for IEAD plot in eV (None if not used)
EEAD_max_energy = 100                        # max energy for IEAD plot in eV (None if not used)

# electric field XT plots
plot_XT_Efield = False  # plot XT plot of z vertical electric field component
XT_colorbar = False    # show color-bars in XT-plots
XT_lower_half = False  # only show bottom half of XT plots

# probe potential over time
plot_local_potential_over_time = True
potential_over_time_locations = [(40, 20), (40, 45), (40, 67)]
potential_over_time_labels = ["Bottom", "Bulk", "Top"]

# root directory
dir_root = os.path.abspath("D:\\UIGEL5_D_Florian\\Voltage_Waveform_Tailoring\\HPEM\\ArCF4O2\\DarkSpaceGeometry"
                               "\\10MHz_Bottom\\10MHz_ConstVoltage_300_1500r1")
# dir_root = os.path.abspath("F:\\Archive\\Voltage_Waveform_Tailoring\\HPEM\\ArCF4O2\\500W_top\\500_2000r2")

# dir_root = os.path.abspath("D:\\UIGEL5_D_Florian\\Voltage_Waveform_Tailoring\\HPEM\\ArCF4O2\\DarkSpaceGeometry\\5MHz_Bottom\\5MHz_ConstVoltage_300_1500res1")

# create config object from locals()
config = Config(locals())
config.print_config()

# ╔═══════════════════════════════════════════════════════════════════════════════════╗ #
# ║   script body                                                                     ║ #
# ╚═══════════════════════════════════════════════════════════════════════════════════╝ #

# initialize case container object
cases = Cases()
# generate case objects by scanning the root folder
cases.scan_directory_for_cases(dir_root)
cases.determine_sweep()
cases.print_info()

# ╔═══════════════════════════════════════════════════════════════════════════════════╗ #
# ║   plotting                                                                        ║ #
# ╚═══════════════════════════════════════════════════════════════════════════════════╝ #

# create directories for figures
path_figures, path_figures_XT_Field, path_figures_IEAD, path_figures_time_varying_potentials \
    = plotting.create_directories(dir_root, config)

# set global plot parameters
plotting.set_plot_globals()

# plot voltages as function of phase
# ----------------------------------------------------------------------------
if not cases.constant_phase:
    # plot DC bias over phase
    xy_phase_bias = cases.get_value_pair("cwaveform_phase", "dc_bias", custom_waveform_only=True)
    plotting.plot_dc_bias_vs_phase(xy_phase_bias, path_figures)
    # plot voltage amplitudes over phase
    xy_phase_voltage1 = cases.get_value_pair("cwaveform_phase", "final_voltages[0]", custom_waveform_only=True)
    xy_phase_voltage2 = cases.get_value_pair("cwaveform_phase", "final_voltages[1]", custom_waveform_only=True)
    plotting.plot_voltages_vs_phase(xy_phase_voltage1, xy_phase_voltage2, xy_phase_bias, path_figures)

# plot geometry
# ----------------------------------------------------------------------------
plotting.plot_geometry(cases, path_figures,
                       plot_local_potential_over_time=plot_local_potential_over_time,
                       potential_over_time_locations=potential_over_time_locations)

# load particle energy-angular distributions from pcmc file and plot
# ----------------------------------------------------------------------------
if plot_EADS or plot_EDFs or plot_ADFs:
    cases.read_pcmc_file_all(config)

    # plot particle energy angular distributions
    if plot_EADS:
        for case in cases:
            for s in plot_EADS_species:
                plotting.plot_ead(case, path_figures_IEAD, iead_max_energy=IEAD_max_energy, plot_species=s)
    if plot_EDFs:
        plotting.plot_edf_compare(cases, path_figures, custom_waveform_only=True, species="ION-TOT")

    if not cases.constant_phase:
        plotting.plot_mean_energies_vs_phase(cases, path_figures, species_plot="ION-TOT")
        plotting.plot_mode_energies_vs_phase(cases, path_figures, species_plot="ION-TOT")

# load and plot time varying voltages
# ----------------------------------------------------------------------------
if plot_local_potential_over_time:

    for case in cases:
        case.get_local_potential_over_time(potential_over_time_locations, potential_over_time_labels)
        plotting.plot_voltages_over_time(case, path_figures_time_varying_potentials)

# create XT-plots of E-field
# ----------------------------------------------------------------------------
if plot_XT_Efield:
    print("\ncreating E-field XT plots...")
    for case in cases:
        plotting.movie2xt(case, path_figures_XT_Field, lower_half=XT_lower_half, do_color_bar=XT_colorbar)

# save data to file
# ----------------------------------------------------------------------------
cases.load_sweep_data()
cases.save_data_to_csv(dir_root)
# TODO:
#  fix determine sweep const voltages and power
