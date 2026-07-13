import argparse
import json
import numpy as np
import uproot
import pprint
import matplotlib.pyplot as plt
from analysis.cut_metadata import cut_metadata
from analysis.create_data_list import create_data_list
from analysis.analyze_files import analyze_files


def plot_mda(data_list, plotcard):
    plot_settings = plotcard["plot"]
    x_key = plot_settings["x_key"]
    x_label = plot_settings["x_label"]
    y_key = plot_settings["y_key"]
    y_label = plot_settings["y_label"]
    legend_title = plot_settings["legend_title"]

    # fig, ax = plt.subplots(3, 1, figsize=(88/inch_to_mm, 130/inch_to_mm), sharex=True, gridspec_kw={"height_ratios": [2, 1, 1]})
    fig, ax = plt.subplots(1, 1, figsize=(88/inch_to_mm, 60/inch_to_mm))

    for line_setting in plot_settings["lines"]:
        line_label = line_setting.pop("plot_label")
        data_points = []
        
        for data_point in data_list:
            if all(data_point[k] == v for k, v in line_setting.items()):
                # data_point["counts_radionuclide"] = counts
                data_points.append(data_point)

        sort_order = sorted(range(len(data_points)), key=lambda i: data_points[i][x_key])
        x_values = [data_points[i][x_key] for i in sort_order]
        y_values = [data_points[i][y_key][0] for i in sort_order]
        y_values_unc = [data_points[i][y_key][1] for i in sort_order]

        ax.errorbar(x_values, y_values, yerr=y_values_unc, fmt=".", ls="-", 
                       markersize=4, lw=1, capsize=2, capthick=1, 
                       label=line_label)

    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    
    ax.legend(frameon=False, fontsize=8, title=legend_title, title_fontsize=8)

    plt.show()


# Parser for adding arguments
parser = argparse.ArgumentParser(prog="plot_mda",
                                 description="Analyze simulation data and plot the MDA",
                                 epilog="Elias Arnqvist, 2026, Uppsala University",
                                 add_help=True)
parser.add_argument("-p", "--plotcard", type=str, required=False, default="plotcards/plotcard_test.json", help="File specifying analysis and plotting settings")
parser.add_argument("-m", "--metadata", type=str, required=False, default="../geant4/output/metadata.json", help="Metadata file")
parser.add_argument("-d", "--data", type=str, required=False, default="../geant4/output/data_pelle", help="Directory of data")
args = parser.parse_args()

# Open the plotcard with analysis and plot settings
plotcard_filepath = args.plotcard
with open(plotcard_filepath, "r") as f:
    plotcard = json.load(f)

# Open the metadata file, describing the properties of the simulation runs
metadata_filepath = args.metadata
with open(metadata_filepath, "r") as f:
    metadata = json.load(f)

# Prepare data filepath
data_path = args.data

metadata = cut_metadata(plotcard, metadata)

data_list = create_data_list(plotcard, metadata)

data_list = analyze_files(data_list, plotcard, metadata, data_path)

pprint.pp(data_list)

# Plot the MDA
plt.close('all')
inch_to_mm = 25.4
# plot_mda(data_list, plotcard)

