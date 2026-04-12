import argparse
import matplotlib.pyplot as plt
import json
import numpy as np


def plot_singles(metadata_backinfo, datacut):
    fig, ax = plt.subplots(1, 1, figsize=(88/inch_to_mm, 60/inch_to_mm))

    back = {}

    for i, (key, value) in enumerate(metadata_backinfo.items()):

        Z = 0
        A = 0

        # Make sure we only plot the data specified in datacut
        cut_included = True
        for j, (key2, value2) in enumerate(datacut.items()):
            if key2 == "plot":
                continue
            elif key2 == "type":
                this_check = (value[key2] == "external_background")
                cut_included *= this_check
            elif key2 == "Z":
                Z = value2
            elif key2 == "A":
                A = value2
            else:
                if value2 == "any":
                    this_check = True
                else:
                    this_check = (value["properties"][key2] == value2)
                cut_included *= this_check

        if not cut_included:
            continue

        # Look at the peak data
        backdata = value["backdata"]["Z" + str(Z) + "_A" + str(A)]["singles"]

        # What to plot on the x-axis
        plot_on_x = datacut["plot"]["x_plot"]        

        for k, E_gamma in enumerate(backdata["E_gamma"]):
            # Data values for this energy
            x_val = value["properties"][plot_on_x]
            counts = backdata["counts_a"][k] + backdata["counts_b"][k]
            pseudo_time = value["properties"]["SURE_pseudo_time"]
            measurement_time = 1*24*60*60
            events = value["properties"]["events"]
            b = counts / pseudo_time
            B = b * measurement_time
            LD = 2.71 + 4.65*np.sqrt(B)
            y_val = LD
            # Use Binomial statistics here instead of Poisson
            db = np.sqrt(counts * (1 - counts/events)) / pseudo_time
            dB = db * measurement_time
            dLD = 4.65 * (0.5/np.sqrt(B)) * dB
            dy_val = dLD

            if E_gamma not in back.keys():
                entry = {"x": [x_val], "y": [y_val], "dy": [dy_val]}
                back[E_gamma] = entry
            else:
                back[E_gamma]["x"].append(x_val)
                back[E_gamma]["y"].append(y_val)
                back[E_gamma]["dy"].append(dy_val)

    # Plot the efficiency
    for l, (key, value) in enumerate(back.items()):
        E_gamma = key
        x = value["x"]
        y = value["y"]
        dy = value["dy"]

        idx = np.argsort(x)
        x_sorted = np.array(x)[idx]
        y_sorted = np.array(y)[idx]
        dy_sorted = np.array(dy)[idx]

        # Generally 2 decimals for ENSDF data, but check for your case!
        line_label = f"{E_gamma:.2f}"
        ax.errorbar(x_sorted, y_sorted, yerr=dy_sorted, fmt=".", ls="-", 
                    markersize=4, lw=1, capsize=2, capthick=1, 
                    label=line_label)
    
    x_label = datacut["plot"]["x_label"]
    ax.set_xlabel(x_label)
    y_label = "$L_D$"
    ax.set_ylabel(y_label)
    legend_title = datacut["plot"]["nuclide"] + str(" (keV)")
    ax.legend(frameon=False, fontsize=8, title=legend_title, title_fontsize=8)

    save_name = "figures/" + datacut["plot"]["save_name"] + "_singles"
    plt.tight_layout(pad = 0.2)
    fig.subplots_adjust(hspace=0, wspace=0)
    plt.savefig(save_name + ".jpg", dpi=300)


def plot_coincidences(metadata_backinfo, datacut):
    fig, ax = plt.subplots(1, 1, figsize=(88/inch_to_mm, 60/inch_to_mm))

    back = {}

    for i, (key, value) in enumerate(metadata_backinfo.items()):

        Z = 0
        A = 0

        # Make sure we only plot the data specified in datacut
        cut_included = True
        for j, (key2, value2) in enumerate(datacut.items()):
            if key2 == "plot":
                continue
            elif key2 == "type":
                this_check = (value[key2] == "external_background")
                cut_included *= this_check
            elif key2 == "Z":
                Z = value2
            elif key2 == "A":
                A = value2
            else:
                if value2 == "any":
                    this_check = True
                else:
                    this_check = (value["properties"][key2] == value2)
                cut_included *= this_check

        if not cut_included:
            continue

        # Look at the peak data
        backdata = value["backdata"]["Z" + str(Z) + "_A" + str(A)]["coincidences"]

        # What to plot on the x-axis
        plot_on_x = datacut["plot"]["x_plot"]        

        for k in range(len(backdata["E_gamma1"])):
            # Coincident energies
            E_gamma1 = backdata["E_gamma1"][k]
            E_gamma2 = backdata["E_gamma2"][k]
            # Data values for this energy
            x_val = value["properties"][plot_on_x]
            counts = backdata["counts_a1b2"][k] + backdata["counts_a2b1"][k]
            pseudo_time = value["properties"]["SURE_pseudo_time"]
            measurement_time = 1*24*60*60
            events = value["properties"]["events"]
            b = counts / pseudo_time
            B = b * measurement_time
            LD = 2.71 + 4.65*np.sqrt(B)
            y_val = LD
            # Use Binomial statistics here instead of Poisson
            db = np.sqrt(counts * (1 - counts/events)) / pseudo_time
            dB = db * measurement_time
            dLD = 4.65 * (0.5/np.sqrt(B)) * dB
            dy_val = dLD

            if k not in back.keys():
                entry = {"x": [x_val], "y": [y_val], "dy": [dy_val], "E_gamma1": E_gamma1, "E_gamma2": E_gamma2}
                back[k] = entry
            else:
                back[k]["x"].append(x_val)
                back[k]["y"].append(y_val)
                back[k]["dy"].append(dy_val)

    # Plot the efficiency
    for l, (key, value) in enumerate(back.items()):
        E_gamma1 = value["E_gamma1"]
        E_gamma2 = value["E_gamma2"]
        x = value["x"]
        y = value["y"]
        dy = value["dy"]

        idx = np.argsort(x)
        x_sorted = np.array(x)[idx]
        y_sorted = np.array(y)[idx]
        dy_sorted = np.array(dy)[idx]

        # Generally 2 decimals for ENSDF data, but check for your case!
        line_label = f"{E_gamma1:.2f}, {E_gamma2:.2f}"
        ax.errorbar(x_sorted, y_sorted, yerr=dy_sorted, fmt=".", ls="-", 
                    markersize=4, lw=1, capsize=2, capthick=1, 
                    label=line_label)
    
    x_label = datacut["plot"]["x_label"]
    ax.set_xlabel(x_label)
    y_label = "$L_D$"
    ax.set_ylabel(y_label)
    legend_title = datacut["plot"]["nuclide"] + str(" (keV, keV)")
    ax.legend(frameon=False, fontsize=8, title=legend_title, title_fontsize=8)

    save_name = "figures/" + datacut["plot"]["save_name"] + "_coincidences"
    plt.tight_layout(pad = 0.2)
    fig.subplots_adjust(hspace=0, wspace=0)
    plt.savefig(save_name + ".jpg", dpi=300)


# Parser for adding arguments
parser = argparse.ArgumentParser(prog="plot_back",
                                 description="Plot the background count rate",
                                 epilog="Elias Arnqvist, 2026, Uppsala University",
                                 add_help=True)
parser.add_argument("-dc", "--datacut", type=str, required=True, default="datacut.json", help="File for only selecting part of the dataset")
parser.add_argument("-mb", "--metadata_backinfo", type=str, required=False, default="metadata_backinfo.json", help="Path to metadata and backinfo file")
args = parser.parse_args()

# Get the metadata and peak counts
metadata_metadata_backinfo_filepath = args.metadata_backinfo
with open(metadata_metadata_backinfo_filepath, "r") as f:
    metadata_backinfo = json.load(f)

# Get the information of how to select what to plot
datacut_filepath = args.datacut
with open(datacut_filepath, "r") as f:
    datacut = json.load(f)

# Plot the efficiency and intensity
plt.close('all')
inch_to_mm = 25.4

plot_singles(metadata_backinfo, datacut)
plot_coincidences(metadata_backinfo, datacut)
