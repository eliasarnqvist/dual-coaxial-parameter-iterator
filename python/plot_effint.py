import argparse
import matplotlib.pyplot as plt
import json
import numpy as np


def plot_singles(metadata_peakinfo, datacut):
    fig, ax = plt.subplots(1, 1, figsize=(88/inch_to_mm, 60/inch_to_mm))

    effint = {}

    for i, (key, value) in enumerate(metadata_peakinfo.items()):

        # Make sure we only plot the data specified in datacut
        cut_included = True
        for j, (key2, value2) in enumerate(datacut.items()):
            if key2 == "plot":
                continue
            elif key2 == "type":
                this_check = (value[key2] == value2)
                cut_included *= this_check
            else:
                if value2 == "any":
                    this_check = True
                else:
                    this_check = (value["properties"][key2] == value2)
                cut_included *= this_check
        
        if not cut_included:
            continue

        # Look at the peak data
        peakdata = value["peakdata"]["singles"]

        # What to plot on the x-axis
        plot_on_x = datacut["plot"]["x_plot"]

        for k, E_gamma in enumerate(peakdata["E_gamma"]):
            # Data values for this energy
            x_val = value["properties"][plot_on_x]
            counts = peakdata["counts_a"][k] + peakdata["counts_b"][k]
            events = value["properties"]["events"]
            y_val = counts / events
            # Use Binomial statistics here instead of Poisson
            dy_val = np.sqrt(counts * (1 - counts/events)) / events

            if E_gamma not in effint.keys():
                entry = {"x": [x_val], "y": [y_val], "dy": [dy_val]}
                effint[E_gamma] = entry
            else:
                effint[E_gamma]["x"].append(x_val)
                effint[E_gamma]["y"].append(y_val)
                effint[E_gamma]["dy"].append(dy_val)

    # Plot the efficiency
    for l, (key, value) in enumerate(effint.items()):
        E_gamma = key
        x = value["x"]
        y = value["y"]
        dy = value["dy"]

        idx = np.argsort(x)
        x_sorted = np.array(x)[idx]
        y_sorted = np.array(y)[idx] * 1e2
        dy_sorted = np.array(dy)[idx] * 1e2

        # Generally 2 decimals for ENSDF data, but check for your case!
        line_label = f"{E_gamma:.2f}"
        ax.errorbar(x_sorted, y_sorted, yerr=dy_sorted, fmt=".", ls="-", 
                    markersize=4, lw=1, capsize=2, capthick=1, 
                    label=line_label)
    
    x_label = datacut["plot"]["x_label"]
    ax.set_xlabel(x_label)
    y_label = r"$I_\gamma \, \varepsilon_\gamma$ (%)"
    ax.set_ylabel(y_label)
    legend_title = datacut["plot"]["nuclide"] + str(" (keV)")
    ax.legend(frameon=False, fontsize=8, title=legend_title, title_fontsize=8)

    save_name = "figures/" + datacut["plot"]["save_name"] + "_singles"
    plt.tight_layout(pad = 0.2)
    fig.subplots_adjust(hspace=0, wspace=0)
    plt.savefig(save_name + ".jpg", dpi=300)


def plot_coincidences(metadata_peakinfo, datacut):
    fig, ax = plt.subplots(1, 1, figsize=(88/inch_to_mm, 60/inch_to_mm))

    effint = {}

    for i, (key, value) in enumerate(metadata_peakinfo.items()):

        # Make sure we only plot the data specified in datacut
        cut_included = True
        for j, (key2, value2) in enumerate(datacut.items()):
            if key2 == "plot":
                continue
            elif key2 == "type":
                this_check = (value[key2] == value2)
                cut_included *= this_check
            else:
                if value2 == "any":
                    this_check = True
                else:
                    this_check = (value["properties"][key2] == value2)
                cut_included *= this_check
        
        if not cut_included:
            continue

        # Look at the peak data
        peakdata = value["peakdata"]["coincidences"]

        # What to plot on the x-axis
        plot_on_x = datacut["plot"]["x_plot"]

        for k in range(len(peakdata["E_gamma1"])):
            # Data values for this energy
            x_val = value["properties"][plot_on_x]
            # Coincident energies
            E_gamma1 = peakdata["E_gamma1"][k]
            E_gamma2 = peakdata["E_gamma2"][k]
            counts = peakdata["counts_a1b2"][k] + peakdata["counts_a2b1"][k]
            events = value["properties"]["events"]
            y_val = counts / events
            # Use Binomial statistics here instead of Poisson
            dy_val = np.sqrt(counts * (1 - counts/events)) / events

            if k not in effint.keys():
                entry = {"x": [x_val], "y": [y_val], "dy": [dy_val], "E_gamma1": E_gamma1, "E_gamma2": E_gamma2}
                effint[k] = entry
            else:
                effint[k]["x"].append(x_val)
                effint[k]["y"].append(y_val)
                effint[k]["dy"].append(dy_val)

    # Plot the efficiency
    for key, value in effint.items():
        E_gamma1 = value["E_gamma1"]
        E_gamma2 = value["E_gamma2"]
        x = value["x"]
        y = value["y"]
        dy = value["dy"]

        idx = np.argsort(x)
        x_sorted = np.array(x)[idx]
        y_sorted = np.array(y)[idx] * 1e2
        dy_sorted = np.array(dy)[idx] * 1e2

        # Generally 2 decimals for ENSDF data, but check for your case!
        line_label = f"{E_gamma1:.2f}, {E_gamma2:.2f}"
        ax.errorbar(x_sorted, y_sorted, yerr=dy_sorted, fmt=".", ls="-", 
                    markersize=4, lw=1, capsize=2, capthick=1, 
                    label=line_label)
    
    x_label = datacut["plot"]["x_label"]
    ax.set_xlabel(x_label)
    y_label = r"$I_{\gamma \gamma} \, \varepsilon_{\gamma \gamma}$ (%)"
    ax.set_ylabel(y_label)
    legend_title = datacut["plot"]["nuclide"] + str(" (keV, keV)")
    ax.legend(frameon=False, fontsize=8, title=legend_title, title_fontsize=8)

    save_name = "figures/" + datacut["plot"]["save_name"] + "_coincidences"
    plt.tight_layout(pad = 0.2)
    fig.subplots_adjust(hspace=0, wspace=0)
    plt.savefig(save_name + ".jpg", dpi=300)


# Parser for adding arguments
parser = argparse.ArgumentParser(prog="plot_effint",
                                 description="Plot the product of efficiency and intensity",
                                 epilog="Elias Arnqvist, 2026, Uppsala University",
                                 add_help=True)
parser.add_argument("-dc", "--datacut", type=str, required=True, default="datacut.json", help="File for only selecting part of the dataset")
parser.add_argument("-mp", "--metadata_peakinfo", type=str, required=False, default="metadata_peakinfo.json", help="Path to metadata and peakinfo file")
args = parser.parse_args()

# Get the metadata and peak counts
metadata_peakinfo_filepath = args.metadata_peakinfo
with open(metadata_peakinfo_filepath, "r") as f:
    metadata_peakinfo = json.load(f)

# Get the information of how to select what to plot
datacut_filepath = args.datacut
with open(datacut_filepath, "r") as f:
    datacut = json.load(f)

# Plot the efficiency and intensity
plt.close('all')
inch_to_mm = 25.4

plot_singles(metadata_peakinfo, datacut)
plot_coincidences(metadata_peakinfo, datacut)
