import argparse
import matplotlib.pyplot as plt
import json
import numpy as np


def plot_mda(metadata_peakinfo, metadata_backinfo, datacut):
    fig, ax = plt.subplots(3, 1, figsize=(88/inch_to_mm, 130/inch_to_mm), sharex=True, gridspec_kw={"height_ratios": [2, 1, 1]})

    # Measurement time in the MDA equation
    # TODO make measurement time a function of half life to get optimal MDA
    measurement_time = 1*24*60*60

    # Dictionary to store efiiciency and backround counts
    effint_back = []
    mda = {}

    # Open both the peak info and back info, they have similar structure
    for i_dict, info_dictionary in enumerate([metadata_peakinfo, metadata_backinfo]):
        for i_info, (key, value) in enumerate(info_dictionary.items()):
            Z = 0
            A = 0

            # Make sure we only plot the data specified in datacut
            cut_included = True
            for i_cut, (key2, value2) in enumerate(datacut.items()):
                # Do not have to look at some of the settings here, like plot settings
                if key2 == "plot_settings" or key2 == "gamma_settings":
                    continue
                # Check the type of the data info
                elif key2 == "type":
                    if info_dictionary == metadata_peakinfo:
                        this_check = (value[key2] == "radionuclides")
                    elif info_dictionary == metadata_backinfo:
                        this_check = (value[key2] == "external_background")
                    else:
                        this_check = False
                    cut_included *= this_check
                # Store the Z and A from the data cut, only for back info as the peak info is handeled below
                elif key2 == "Z" and info_dictionary == metadata_backinfo:
                    Z = value2
                elif key2 == "A" and info_dictionary == metadata_backinfo:
                    A = value2
                # Check the rest of the data cut compared with the info files
                else:
                    if value2 == "any":
                        this_check = True
                    # Special case only for the background info
                    elif info_dictionary == metadata_peakinfo and key2 == "background_file":
                        this_check = True
                    # General check of cut and properties of the run
                    else:
                        this_check = (value["properties"][key2] == value2)
                    cut_included *= this_check

            # If the include cut is not true we continue withy the next entry of the peak or background info file
            if not cut_included:
                continue
            
            # Get the relevant singles data
            for E_gamma in datacut["gamma_settings"]["singles"]["E_gamma"]:
                if info_dictionary == metadata_peakinfo:
                    peakdata = value["peakdata"]["singles"]

                    index = peakdata["E_gamma"].index(E_gamma)
                    counts = peakdata["counts_a"][index] + peakdata["counts_b"][index]
                    events = value["properties"]["events"]

                    y_val = counts / events * 1e2 # in units of percent
                    dy_val = np.sqrt(counts * (1 - counts/events)) / events * 1e2 # in untis of percent

                    plot_on_x = datacut["plot_settings"]["x_plot"]
                    x_val = value["properties"][plot_on_x]

                    entry = {"type": "radionuclides",
                             "detector": "singles",
                             "label": f"{E_gamma:.2f}",
                             "x_val": x_val,
                             "effint": y_val,
                             "effint_unc": dy_val}
                    effint_back.append(entry)

                elif info_dictionary == metadata_backinfo:
                    backdata = value["backdata"]["Z" + str(Z) + "_A" + str(A)]["singles"]

                    index = backdata["E_gamma"].index(E_gamma)
                    counts = backdata["counts_a"][index] + backdata["counts_b"][index]
                    events = value["properties"]["events"]

                    pseudo_time = value["properties"]["SURE_pseudo_time"]
                    b = counts / pseudo_time # background count rate
                    B = b * measurement_time # measured background counts
                    LD = 2.71 + 4.65*np.sqrt(B) # detection limit (Currie for now, Stapleton later?)
                    y_val = LD
                    db = np.sqrt(counts * (1 - counts/events)) / pseudo_time
                    dB = db * measurement_time
                    dLD = 4.65 * (0.5/np.sqrt(B)) * dB
                    dy_val = dLD

                    plot_on_x = datacut["plot_settings"]["x_plot"]
                    x_val = value["properties"][plot_on_x]

                    entry = {"type": "background",
                             "detector": "singles",
                             "label": f"{E_gamma:.2f}",
                             "x_val": x_val,
                             "LD": y_val,
                             "LD_unc": dy_val}
                    effint_back.append(entry)

            # Get the relevant coincidences data
            for E_gamma1, E_gamma2 in zip(datacut["gamma_settings"]["coincidences"]["E_gamma1"], datacut["gamma_settings"]["coincidences"]["E_gamma2"]):
                # print(E_gamma1, E_gamma2)

                if info_dictionary == metadata_peakinfo:
                    peakdata = value["peakdata"]["coincidences"]

                    index = np.where((np.array(peakdata["E_gamma1"]) == E_gamma1) & (np.array(peakdata["E_gamma2"]) == E_gamma2))[0][0]
                    counts = peakdata["counts_a1b2"][index] + peakdata["counts_a2b1"][index]
                    events = value["properties"]["events"]

                    y_val = counts / events * 1e2 # in units of percent
                    dy_val = np.sqrt(counts * (1 - counts/events)) / events * 1e2 # in untis of percent

                    plot_on_x = datacut["plot_settings"]["x_plot"]
                    x_val = value["properties"][plot_on_x]

                    # if value["properties"]["select_ntype_instead_of_ptype"] == True:
                    #     label = f"n-type HPGe"
                    # elif value["properties"]["select_ntype_instead_of_ptype"] == False:
                    #     label = f"p-type HPGe"
                    label = f"{E_gamma1:.2f}, {E_gamma2:.2f}"
                    
                    entry = {"type": "radionuclides",
                             "detector": "coincidences",
                             "label": label,
                             "x_val": x_val,
                             "effint": y_val,
                             "effint_unc": dy_val}
                    effint_back.append(entry)
                
                elif info_dictionary == metadata_backinfo:
                    backdata = value["backdata"]["Z" + str(Z) + "_A" + str(A)]["coincidences"]

                    index = np.where((np.array(backdata["E_gamma1"]) == E_gamma1) & (np.array(backdata["E_gamma2"]) == E_gamma2))[0][0]

                    events = value["properties"]["events"]
                    ROI = backdata["ROI"]
                    ROI_standard = 5
                    if ROI != ROI_standard:
                        counts_raw = backdata["counts_a1b2"][index] + backdata["counts_a2b1"][index]
                        counts = counts_raw / ((ROI/ROI_standard)**2)
                        counts_raw_unc = np.sqrt(counts_raw * (1 - counts_raw/events))
                        counts_unc = counts_raw_unc / ((ROI/ROI_standard)**2)
                    else:
                        counts = backdata["counts_a1b2"][index] + backdata["counts_a2b1"][index]
                        counts_unc = np.sqrt(counts * (1 - counts/events))
                    pseudo_time = value["properties"]["SURE_pseudo_time"]
                    b = counts / pseudo_time # background count rate
                    B = b * measurement_time # measured background counts
                    LD = 2.71 + 4.65*np.sqrt(B) # detection limit (Currie for now, Stapleton later?)
                    y_val = LD
                    db = counts_unc / pseudo_time
                    dB = db * measurement_time
                    dLD = 4.65 * (0.5/np.sqrt(B)) * dB
                    dy_val = dLD

                    plot_on_x = datacut["plot_settings"]["x_plot"]
                    x_val = value["properties"][plot_on_x]

                    # if value["properties"]["select_ntype_instead_of_ptype"] == True:
                    #     label = f"n-type HPGe"
                    # elif value["properties"]["select_ntype_instead_of_ptype"] == False:
                    #     label = f"p-type HPGe"
                    label = f"{E_gamma1:.2f}, {E_gamma2:.2f}"

                    entry = {"type": "background",
                             "detector": "coincidences",
                             "label": label,
                             "x_val": x_val,
                             "LD": y_val,
                             "LD_unc": dy_val}
                    effint_back.append(entry)
            
    # print(effint_back)

    # Sort all the values based on the label which will be plotted later
    for effint_back_dict in effint_back:
        label = effint_back_dict["label"]
        if label not in mda.keys():
            mda[label] = {"effint_x": np.array([]), "effint": np.array([]), "effint_unc": np.array([]),
                          "LD_x": np.array([]), "LD": np.array([]), "LD_unc": np.array([]),
                          "detector": effint_back_dict["detector"]}

        if effint_back_dict["type"] == "radionuclides":
            mda[label]["effint_x"] = np.append(mda[label]["effint_x"], effint_back_dict["x_val"])
            mda[label]["effint"] = np.append(mda[label]["effint"], effint_back_dict["effint"])
            mda[label]["effint_unc"] = np.append(mda[label]["effint_unc"], effint_back_dict["effint_unc"])

        elif effint_back_dict["type"] == "background":
            mda[label]["LD_x"] = np.append(mda[label]["LD_x"], effint_back_dict["x_val"])
            mda[label]["LD"] = np.append(mda[label]["LD"], effint_back_dict["LD"])
            mda[label]["LD_unc"] = np.append(mda[label]["LD_unc"], effint_back_dict["LD_unc"])

    # Sort the values based on the x-value specified and caluclate the MDA
    for key, value in mda.items():
        effint_i_sorted = np.argsort(value["effint_x"])

        mda[key]["effint_x"] = value["effint_x"][effint_i_sorted]
        mda[key]["effint"] = value["effint"][effint_i_sorted]
        mda[key]["effint_unc"] = value["effint_unc"][effint_i_sorted]

        back_i_sorted = np.argsort(value["LD_x"])

        mda[key]["LD_x"] = value["LD_x"][back_i_sorted]
        mda[key]["LD"] = value["LD"][back_i_sorted]
        mda[key]["LD_unc"] = value["LD_unc"][back_i_sorted]

        effint = mda[key]["effint"]
        effint_unc = mda[key]["effint_unc"]
        LD = mda[key]["LD"]
        LD_unc = mda[key]["LD_unc"]
        MDA = LD / (effint * measurement_time)
        MDA_unc = np.sqrt((1/(effint*measurement_time) * LD_unc)**2 + (LD/(effint**2*measurement_time) * effint_unc)**2)

        mda[key]["MDA_x"] = mda[key]["effint_x"] # Really needed?
        mda[key]["MDA"] = MDA
        mda[key]["MDA_unc"] = MDA_unc

    # Plot the effint, back, and mda
    for key, value in mda.items():
        label = key
        x = value["MDA_x"]

        effint = value["effint"]
        effint_unc = value["effint_unc"]

        LD = value["LD"]
        LD_unc = value["LD_unc"]

        MDA = value["MDA"]
        MDA_unc = value["MDA_unc"]

        # Plot the 
        ax[0].errorbar(x, MDA, yerr=MDA_unc, fmt=".", ls="-", 
                        markersize=4, lw=1, capsize=2, capthick=1, 
                        label=label)
        
        ax[1].errorbar(x, effint, yerr=effint_unc, fmt=".", ls="-", 
                        markersize=4, lw=1, capsize=2, capthick=1)

        ax[2].errorbar(x, LD, yerr=LD_unc, fmt=".", ls="-", 
                        markersize=4, lw=1, capsize=2, capthick=1)
    
    ax[0].set_zorder(3)
    ax[1].set_zorder(2)
    ax[2].set_zorder(1)

    ax[0].set_ylabel("MDA (Bq)")
    ax[1].set_ylabel(r"$I_{\gamma \gamma} \, \varepsilon_{\gamma \gamma}$, $I_{\gamma} \, \varepsilon_{\gamma}$ (%)")
    ax[2].set_ylabel(r"$L_D$")
    x_label = datacut["plot_settings"]["x_label"]
    ax[2].set_xlabel(x_label)

    legend_title = datacut["plot_settings"]["nuclide"]
    ax[0].legend(frameon=False, fontsize=8, title=legend_title, title_fontsize=8)

    ax[0].set_yscale("log")
    ax[1].set_yscale("log")
    ax[2].set_yscale("log")

    save_name = "figures/" + datacut["plot_settings"]["save_name"]
    plt.tight_layout(pad = 0.2)
    fig.subplots_adjust(hspace=0, wspace=0)
    plt.savefig(save_name + ".jpg", dpi=300)

    plt.show()


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

plot_mda(metadata_peakinfo, metadata_backinfo, datacut)
