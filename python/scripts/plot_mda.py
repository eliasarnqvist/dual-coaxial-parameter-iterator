import argparse
import json
import time
import subprocess
import uuid
import os
import numpy as np
import uproot
from numba import jit


# step 1: cut the metadata
# OK!

def cut_metadata(plotcard, metadata):
    # Use the selection cuts in the plotcard to remove irrelevant metadata
    for key in list(metadata):
        value = metadata[key]
        
        keep_this_key = True

        # Check that the geometry matches (for all runs)
        for geometry_setting in ["detector_type", "detector_diameter", "detector_length", "detector_source_distance", "source_type"]:
            if len(plotcard["geometry"][geometry_setting]) == 1 and plotcard["geometry"][geometry_setting][0] == "any":
                pass
            elif value["properties"][geometry_setting] not in plotcard["geometry"][geometry_setting]:
                keep_this_key = False
            else:
                pass

        # Check that the radionuclide/background settings match (only for applicable runs)
        if value["type"] == "radionuclides":
            ZA = [value["properties"]["Z"], value["properties"]["A"]]
            if ZA not in plotcard["radionuclides"]["ZAs"]:
                keep_this_key = False
        elif value["type"] == "background":
            if value["properties"]["background_file"] not in plotcard["background"]["background_file"]:
                keep_this_key = False
        elif value["type"] == "filter":
            ZA = [value["properties"]["Z"], value["properties"]["A"]]
            ZAs_filter = [plotcard["filter"]["ZAs_Bq_filter"][i][0:2] for i in range(len(plotcard["filter"]["ZAs_Bq_filter"]))]
            ZAs_exposure = [plotcard["filter"]["ZAs_Bq_filter"][i][0:2] for i in range(len(plotcard["filter"]["ZAs_Bq_filter"]))]
            if (ZA not in ZAs_filter) and (ZA not in ZAs_exposure):
                keep_this_key = False
        else:
            print("BAD RUN TYPE!")
            raise ValueError
        
        # Remove the key from the metadata dictionary if it is not interesting
        if keep_this_key == False:
            metadata.pop(key)
    
    # Return the shortened metadata dictionary
    return metadata


# step 2: open and analyze the data: ROI, effint, back, filter_back (both intrisic and exposure!), and MDA
# analyze: input: file, output effint

def analyze_file(plotcard, metadata, data_path):
    print("Starting analysis")

    # Need to iterate through the relevant files
    for i, (key, value) in enumerate(metadata.items()):
        run_type = value["type"]

        data_filename = value["filename"]
        data_filesize = value["file_size"] / (1024**2)
        print("\tOpening file " + str(i+1) + " out of " + str(len(metadata)) + ": " + str(data_filename) + " (" + run_type + f" / {data_filesize:.2f} mb)...")

        data_filepath = data_path + "/" + data_filename
        data_file = uproot.open(data_filepath)
        # Single gamma rays
        tree_detector_a = data_file["E_a_list"]
        tree_detector_b = data_file["E_b_list"]
        events_single_a = tree_detector_a["energy_a"].array(library="np") * 1e3
        events_single_b = tree_detector_b["energy_b"].array(library="np") * 1e3
        # Coincidence gamma rays
        tree_detector_ab = data_file["E_ab_list"]
        events_coincidence_a = tree_detector_ab["energy_a"].array(library="np") * 1e3
        events_coincidence_b = tree_detector_ab["energy_b"].array(library="np") * 1e3
        # At the moment the data loading above is the bottleneck...

        radionuclide_counts_singles = []
        radionuclide_counts_coincidences = []

        # Need to iterate through every peak of interest, which means first iterating though the radionuclides
        for j, (Z, A) in enumerate(plotcard["radionuclides"]["ZAs"]):

            print("\t\tAnalyzing radionuclide " + str(j+1) + " out of " + str(len(plotcard["radionuclides"]["ZAs"])) + f": Z={Z}, A={A}")

            counts_singles = []
            counts_coincidences = []

            # Need to iterate through the gamma singles peaks
            for k, E_gamma in enumerate(plotcard["radionuclides"]["gammas_singles"][j]):
                # Need to handle the combined peaks later
                if E_gamma == "combined":
                    continue
                else:
                    pass
                
                ROI_polynial_coeffs = plotcard["analysis"]["ROI_width_polynomial"]
                ROI_width = np.polyval(ROI_polynial_coeffs, E_gamma)

                counts_a = ROI_analysis_1D(events_single_a, E_gamma, ROI_width)
                counts_b = ROI_analysis_1D(events_single_b, E_gamma, ROI_width)
                print("\t\t\tGamma ray: " + str(E_gamma) + " keV, counts a: " + str(counts_a) + ", counts b: " + str(counts_b))

                counts_singles.append([counts_a, counts_b])

            # Need to iterate through the gamma coincidences peaks
            for k, E_gammas in enumerate(plotcard["radionuclides"]["gammas_coincidences"][j]):
                if E_gammas == "combined":
                    continue
                else:
                    E_gamma1, E_gamma2 = E_gammas
                    pass

                ROI_polynial_coeffs = plotcard["analysis"]["ROI_width_polynomial"]
                ROI_width1 = np.polyval(ROI_polynial_coeffs, E_gamma1)
                ROI_width2 = np.polyval(ROI_polynial_coeffs, E_gamma2)

                counts_a1b2 = ROI_analysis_2D(events_coincidence_a, events_coincidence_b, E_gamma1, E_gamma2, ROI_width1, ROI_width2)
                counts_b1a2 = ROI_analysis_2D(events_coincidence_b, events_coincidence_a, E_gamma1, E_gamma2, ROI_width1, ROI_width2)
                print("\t\t\tGamma ray 1: " + str(E_gamma1) + " keV, gamma ray 2: " + str(E_gamma2) + " keV , counts a1b2: " + str(counts_a1b2) + ", counts b2a1: " + str(counts_b1a2))

                counts_coincidences.append([counts_a1b2, counts_b1a2])
            
            radionuclide_counts_singles.append(counts_singles)
            radionuclide_counts_coincidences.append(counts_coincidences)
        
        value["counts_singles"] = radionuclide_counts_singles
        value["counts_coincidences"] = radionuclide_counts_coincidences
        plotcard["metadata"][key] = value
        


    # Put the counts and metadata in the plotcard now


    # Do combined now! For MDA!!! Both 1D 2D
            # if "combined" in plotcard["radionuclides"]["gammas_singles"][j]:
            #     print("do combined also")

    # result should be LD, effint, mda, x value, and line label

    # need: gamma energies from plotcard
    # get energy -> get ROI
    # go through data look for counts in ROI: single and coincidence


    return


def ROI_analysis_1D(events, E_gamma, ROI_width):
    # Energy distance is half of the ROI size
    dE = ROI_width/2

    # Count peak counts
    cond = np.logical_and(events > E_gamma-dE, events < E_gamma+dE)
    counts = cond.sum()
    return int(counts)



def ROI_analysis_2D(events_a, events_b, E_gamma1, E_gamma2, ROI_width1, ROI_width2):
    # Energy distance is half of the ROI size
    dE_1 = ROI_width1/2
    dE_2 = ROI_width2/2

    # Count peak counts
    cond_1 = np.logical_and(events_a > E_gamma1-dE_1, events_a < E_gamma1+dE_1)
    cond_2 = np.logical_and(events_b > E_gamma2-dE_2, events_b < E_gamma2+dE_2)
    cond = np.logical_and(cond_1, cond_2)
    counts = cond.sum()
    return int(counts)


def calculate_LD():

    return


def calculate_effint():

    return


# step 3: mda calculation

# step 4: plot and save








# Parser for adding arguments
parser = argparse.ArgumentParser(prog="plot_mda",
                                 description="Analyze simulation data and plot the MDA",
                                 epilog="Elias Arnqvist, 2026, Uppsala University",
                                 add_help=True)
parser.add_argument("-p", "--plotcard", type=str, required=False, default="../plotcards/plotcard_test.json", help="File specifying analysis and plotting settings")
parser.add_argument("-m", "--metadata", type=str, required=False, default="../../geant4/output/metadata.json", help="Metadata file")
parser.add_argument("-d", "--data", type=str, required=False, default="../../geant4/output/data_pelle", help="Directory of data")
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

analyze_file(plotcard, metadata, data_path)


# print(plotcard)

# print(metadata)


