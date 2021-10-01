import matplotlib.pyplot as plt


# plotting presets
def set_plot_globals():
    plt.rcParams['agg.path.chunksize'] = 10000
    #plt.rcParams['mathtext.fontset'] = 'stix'
    #plt.rcParams['font.family'] = 'STIXGeneral'
    plt.rcParams['font.size'] = 16

def plot_dc_bias_over_phase(xy):
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
    plt.show()
