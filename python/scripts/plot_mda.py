import argparse
import json
import time
import subprocess
import uuid
import os
import numpy as np
import uproot
import pprint


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


def analyze_file(plotcard, metadata, data_path):
    print("Starting analysis")

    # Need to iterate through the relevant files
    for i, (key, value) in enumerate(metadata.items()):

        # NOTE just for testing
        if i > 5:
            continue

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

    return plotcard


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


def calculate_quantities(plotcard):

    for i, (key, value) in enumerate(plotcard["metadata"].items()):

        if value["type"] == "radionuclides":
            # If the type is radionuclides, the efficiency and intensity should be calculated

            # Use both detectors as a bigger single detector
            effint_singles = []
            for radionuclide_counts in value["counts_singles"]:
                radionuclide_effint = []
                for gamma_counts in radionuclide_counts:
                    total_counts = sum(gamma_counts)
                    events = value["properties"]["events"]

                    effint, effint_unc = caluclate_effint(total_counts, events)

                    radionuclide_effint.append([effint, effint_unc])
                effint_singles.append(radionuclide_effint)
            
            # Use standard coincidence detection (no addback)
            effint_coincidences = []
            for radionuclide_counts in value["counts_coincidences"]:
                radionuclide_effint = []
                for gamma_counts in radionuclide_counts:
                    total_counts = sum(gamma_counts)
                    events = value["properties"]["events"]

                    effint, effint_unc = caluclate_effint(total_counts, events)

                    radionuclide_effint.append([effint, effint_unc])
                effint_coincidences.append(radionuclide_effint)
            
            # Add the efficiencies to the metadata
            value["effint_singles"] = effint_singles
            value["effint_coincidences"] = effint_coincidences

        elif value["type"] == "background":
            # If the type is background, the detection limit should be calculated

            # Use both detectors as a bigger single detector
            LD_singles = []
            for radionuclide_counts in value["counts_singles"]:
                radionuclide_LD = []
                for gamma_counts in radionuclide_counts:
                    total_counts = sum(gamma_counts)
                    events = value["properties"]["events"]
                    pseudo_time = value["properties"]["SURE_pseudo_time"]
                    measurement_time = plotcard["analysis"]["measurement_time_hours"]

                    LD, LD_unc = calculate_LD(total_counts, events, pseudo_time, measurement_time)

                    radionuclide_LD.append([LD, LD_unc])
                LD_singles.append(radionuclide_LD)

            # Use standard coincidence detection (no addback)
            LD_coincidences = []
            for radionuclide_counts in value["counts_coincidences"]:
                radionuclide_LD = []
                for gamma_counts in radionuclide_counts:
                    total_counts = sum(gamma_counts)
                    events = value["properties"]["events"]
                    pseudo_time = value["properties"]["SURE_pseudo_time"]
                    measurement_time = plotcard["analysis"]["measurement_time_hours"]

                    LD, LD_unc = calculate_LD(total_counts, events, pseudo_time, measurement_time)

                    radionuclide_LD.append([LD, LD_unc])
                LD_coincidences.append(radionuclide_LD)


            # Add the efficiencies to the metadata
            value["LD_singles"] = LD_singles
            value["LD_coincidences"] = LD_coincidences

        elif value["type"] == "filter":
            # Pass for now, TODO later
            pass
        else:
            pass

        # Update the metadata in the plotcad
        plotcard["metadata"][key] = value

    # pprint.pp(plotcard)

    return plotcard


def caluclate_effint(total_counts, events):
    # Calculate efficiency and intensity
    effint = total_counts / events

    # Use binomial instead of poisson approximation
    effint_unc = np.sqrt(total_counts * (1 - total_counts/events)) / events

    return effint, effint_unc


# TODO implement proper detection limit formula that works for low counts
# TODO implement uncertainty of B when B is zero (one-sided uncertainty)
def calculate_LD(total_counts, events, pseudo_time, measurement_time):
    # Background count rate
    b = total_counts / pseudo_time
    # Background counts during measurement time
    B = b * (measurement_time * 3600)
    if B == 0:
        print("Encountered 0 background!!!")
        raise ValueError
    # Detection limit according to Currie
    LD = 2.71 + 4.65*np.sqrt(B)

    # Calculate the uncertainty
    b_unc = np.sqrt(total_counts * (1 - total_counts/events)) / pseudo_time
    B_unc = b_unc * (measurement_time * 3600)
    LD_unc = 4.65 * (0.5/np.sqrt(B)) * B_unc

    return LD, LD_unc
            

def plot_mda():
    # Now the MDA can be caluclated in the next step
    return


def calculate_mda(LD, effint, tM):
    # Calculate the minimum detectable activity
    mda = LD / (effint * tM)

    return mda


def calculate_combined_mda():

    # Check if the combined MDA should be calculated
    # if "combined" in plotcard["radionuclides"]["gammas_singles"][j]:
    #     print("do combined also")
    # This should be for every radionucldie... need to iterate through ZAs
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

plotcard = analyze_file(plotcard, metadata, data_path)

plotcard = calculate_quantities(plotcard)



# print(plotcard)

# print(metadata)


