import csv
import matplotlib.pyplot as plt
import plotting

# set global plot parameters
plotting.set_plot_globals()

root = "D:\\UIGEL5_D_Florian\\Voltage_Waveform_Tailoring\\Figures\\Low_Frequency_Limit\\const_volt_mean energy comparison\\"

paths = [
    root+"1MHzdata.csv",
    root+"2MHzdata.csv",
    root+"5MHzdata.csv",
    root+"10MHzdata.csv"]

labels = [
    "1MHz",
    "2MHz",
    "5MHz",
    "10MHz"]

labels_posx = [0, 45, 90, 135]


name_var1 = "phase"
name_var2 = "mean energy"

plotting.plot_compare_sweeps_mean_energy(paths, name_var1, name_var2, labels, labels_posx)
