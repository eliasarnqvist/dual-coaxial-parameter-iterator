import argparse
import matplotlib.pyplot as plt
import json
import numpy as np


def plot_mda(metadata_peakinfo, metadata_backinfo, datacut, plot_coincidences=True):
    fig, ax = plt.subplots(3, 1, figsize=(88/inch_to_mm, 140/inch_to_mm), sharex=True, gridspec_kw={"height_ratios": [2, 1, 1]})

    measurement_time = 1*24*60*60

    effint_back = {}

    for info_dictionary in [metadata_peakinfo, metadata_backinfo]:
        for i, (key, value) in enumerate(info_dictionary.items()):
            Z = 0
            A = 0

            # Make sure we only plot the data specified in datacut
            cut_included = True
            for j, (key2, value2) in enumerate(datacut.items()):
                if key2 == "plot":
                    continue
                elif key2 == "type":
                    if info_dictionary == metadata_peakinfo:
                        this_check = (value[key2] == "radionuclides")
                    elif info_dictionary == metadata_backinfo:
                        this_check = (value[key2] == "external_background")
                    else:
                        this_check = False
                    cut_included *= this_check
                elif key2 == "Z" and info_dictionary == metadata_backinfo:
                    Z = value2
                elif key2 == "A" and info_dictionary == metadata_backinfo:
                    A = value2
                else:
                    if value2 == "any" or (key2 == "background_file" and info_dictionary == metadata_peakinfo):
                        this_check = True
                    else:
                        this_check = (value["properties"][key2] == value2)
                    cut_included *= this_check

            if not cut_included:
                continue

            # Look at peak data and back data
            if plot_coincidences:
                gamma_category = "coincidences"
            else:
                gamma_category = "singles"
            
            if info_dictionary == metadata_peakinfo:
                peakdata = value["peakdata"][gamma_category]
                gamma_iterator = range(len(peakdata["E_gamma1"]))
            elif info_dictionary == metadata_backinfo:
                backdata = value["backdata"]["Z" + str(Z) + "_A" + str(A)][gamma_category]
                gamma_iterator = range(len(backdata["E_gamma1"]))

            # What to plot on the x-axis
            plot_on_x = datacut["plot"]["x_plot"]

            for k in gamma_iterator:
                if info_dictionary == metadata_peakinfo:
                    # Coincident energies
                    E_gamma1 = peakdata["E_gamma1"][k]
                    E_gamma2 = peakdata["E_gamma2"][k]
                elif info_dictionary == metadata_backinfo:
                    # Coincident energies
                    E_gamma1 = backdata["E_gamma1"][k]
                    E_gamma2 = backdata["E_gamma2"][k]

                # Data values for this energy
                x_val = value["properties"][plot_on_x]

                key_name = (E_gamma1, E_gamma2)
                entry = effint_back.setdefault(key_name, {"x_effint": [], "effint": [], "deffint": [], "x_ld": [], "ld": [], "dld": []})

                if info_dictionary == metadata_peakinfo:
                    # Effint
                    counts = peakdata["counts_a1b2"][k] + peakdata["counts_a2b1"][k]
                    events = value["properties"]["events"]
                    y_val = counts / events * 1e2 # in units of percent
                    dy_val = np.sqrt(counts * (1 - counts/events)) / events * 1e2 # in untis of percent

                    entry["x_effint"].append(x_val)
                    entry["effint"].append(y_val)
                    entry["deffint"].append(dy_val)
                elif info_dictionary == metadata_backinfo:                    
                    # Back
                    counts = backdata["counts_a1b2"][k] + backdata["counts_a2b1"][k]
                    events = value["properties"]["events"]
                    pseudo_time = value["properties"]["SURE_pseudo_time"]
                    b = counts / pseudo_time # background count rate
                    B = b * measurement_time # background measured counts
                    LD = 2.71 + 4.65*np.sqrt(B) # detection limit (Currie for now, Stapleton later)
                    y_val = LD
                    db = np.sqrt(counts * (1 - counts/events)) / pseudo_time
                    dB = db * measurement_time
                    dLD = 4.65 * (0.5/np.sqrt(B)) * dB
                    dy_val = dLD

                    entry["x_ld"].append(x_val)
                    entry["ld"].append(y_val)
                    entry["dld"].append(dy_val)
    
    # Plot the effint, back, and mda
    for (E_gamma1, E_gamma2), value in effint_back.items():
        x_effint = value["x_effint"]
        effint = value["effint"]
        deffint = value["deffint"]
        idx_effint = np.argsort(x_effint)
        x_effint_sorted = np.array(x_effint)[idx_effint]
        effint_sorted = np.array(effint)[idx_effint]
        deffint_sorted = np.array(deffint)[idx_effint]

        x_ld = value["x_ld"]
        ld = value["ld"]
        dld = value["dld"]
        idx_back = np.argsort(x_ld)
        x_ld_sorted = np.array(x_ld)[idx_back]
        ld_sorted = np.array(ld)[idx_back]
        dld_sorted = np.array(dld)[idx_back]

        mda = ld_sorted / (effint_sorted * measurement_time)
        dmda = np.sqrt((1/(effint_sorted*measurement_time) * dld_sorted)**2 + (ld_sorted/(effint_sorted**2*measurement_time) * deffint_sorted)**2)

        # Generally 2 decimals for ENSDF data, but check for your case!
        line_label = f"{E_gamma1:.2f}, {E_gamma2:.2f}"

        ax[0].errorbar(x_effint_sorted, mda, yerr=dmda, fmt=".", ls="-", 
                        markersize=4, lw=1, capsize=2, capthick=1, 
                        label=line_label)
        
        ax[1].errorbar(x_effint_sorted, effint_sorted, yerr=deffint_sorted, fmt=".", ls="-", 
                        markersize=4, lw=1, capsize=2, capthick=1)

        ax[2].errorbar(x_effint_sorted, ld_sorted, yerr=dld_sorted, fmt=".", ls="-", 
                        markersize=4, lw=1, capsize=2, capthick=1)
    
    ax[0].set_zorder(3)
    ax[1].set_zorder(2)
    ax[2].set_zorder(1)

    ax[0].set_yscale("log")
    # ax[1].set_yscale("log")

    x_label = datacut["plot"]["x_label"]
    ax[2].set_xlabel(x_label)
    ax[0].set_ylabel("MDA (Bq)")
    ax[1].set_ylabel(r"$I_{\gamma \gamma} \, \varepsilon_{\gamma \gamma}$ (%)")
    ax[2].set_ylabel(r"$L_D$")
    legend_title = datacut["plot"]["nuclide"] + str(" (keV, keV)")
    ax[0].legend(frameon=False, fontsize=8, title=legend_title, title_fontsize=8)

    save_name = "figures/" + datacut["plot"]["save_name"] + "_coincidences"
    plt.tight_layout(pad = 0.2)
    fig.subplots_adjust(hspace=0, wspace=0)
    plt.savefig(save_name + ".jpg", dpi=300)

# Parser for adding arguments
parser = argparse.ArgumentParser(prog="plot_mda",
                                 description="Plot the minimum detectable activity",
                                 epilog="Elias Arnqvist, 2026, Uppsala University",
                                 add_help=True)
parser.add_argument("-dc", "--datacut", type=str, required=True, default="datacut.json", help="File for only selecting part of the dataset")
parser.add_argument("-mp", "--metadata_peakinfo", type=str, required=False, default="metadata_peakinfo.json", help="Path to metadata and peakinfo file")
parser.add_argument("-mb", "--metadata_backinfo", type=str, required=False, default="metadata_backinfo.json", help="Path to metadata and backinfo file")
args = parser.parse_args()

# Get the metadata and peak counts
metadata_peakinfo_filepath = args.metadata_peakinfo
with open(metadata_peakinfo_filepath, "r") as f:
    metadata_peakinfo = json.load(f)

# Get the metadata and background counts
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

# plot_singles(metadata_backinfo, datacut)
# plot_coincidences(metadata_backinfo, datacut)

plot_mda(metadata_peakinfo, metadata_backinfo, datacut, plot_coincidences=True)
# plot_mda(metadata_peakinfo, metadata_backinfo, datacut, plot_coincidences=False)

