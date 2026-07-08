import argparse
import json
import time
import subprocess
import uuid
import os
import numpy as np
import uproot
import pprint
# from analysis.cut_metadata import cut_metadata


def cut_metadata(plotcard, metadata):
    # Cut the metadata according to the cuts in the plotcard

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
            raise ValueError
        
        # Remove the key from the metadata dictionary if it is not interesting
        if keep_this_key == False:
            metadata.pop(key)
    
    # Return the shortened metadata dictionary
    return metadata


def create_data_list(plotcard, metadata):
    # Make a "list mode" data storage scheme
    # A list with dictionaries
    # Every list entry is a specific data point with a specific detector design, radionuclide, and gamma ray(s)
    # The dictionary specifies which properties this data point has
    # Calculated quantities (like the MDAs) are added to this list later on

    data_list = []
    
    # First use the metadata to supply the geometry and radionuclide information
    for key, value in metadata.items():
        if value["type"] == "radionuclides":
            data_point = {
                "detector_type": value["properties"]["detector_type"],
                "detector_diameter": value["properties"]["detector_diameter"],
                "detector_length": value["properties"]["detector_length"],
                "detector_source_distance": value["properties"]["detector_source_distance"],
                "source_type": value["properties"]["source_type"],
                "ZA": [value["properties"]["Z"], value["properties"]["A"]],
            }
            if data_point not in data_list:
                data_list.append(data_point)
            else:
                # Somehow there is a duplicate!
                raise ValueError
        else:
            # Not metadata for a radionuclide
            pass
    
    # Next use the plotcard to supply the background file information
    new_data_list = []
    background_files = plotcard["background"]["background_file"]
    for data_point in data_list:
        for background_file in background_files:
            new_data_point = data_point.copy() # important with .copy() here!
            new_data_point["background_file"] = background_file
            new_data_list.append(new_data_point)
    data_list = new_data_list

    # Next use the plotcard to supply the gamma information
    new_data_list = []
    ZAs = plotcard["radionuclides"]["ZAs"]
    for data_point in data_list:
        for i, ZA in enumerate(ZAs):
            if data_point["ZA"] == ZA:
                # Check which energies are relevant for this radionuclide
                gammas_singles = plotcard["radionuclides"]["gammas_singles"][i]
                gammas_coincidences = plotcard["radionuclides"]["gammas_coincidences"][i]

                # Need a new data point for every single gamma ray energy
                for gamma_single in gammas_singles:
                    new_data_point = data_point.copy() # important with .copy() here!
                    new_data_point["analysis_type"] = "singles"
                    new_data_point["Egamma"] = gamma_single
                    new_data_list.append(new_data_point)
                
                # Need a new data point for every pair of coincidence gamma rays
                for gamma_coincidence in gammas_coincidences:
                    new_data_point = data_point.copy() # important with .copy() here!
                    new_data_point["analysis_type"] = "coincidences"
                    new_data_point["Egamma1_Egamma2"] = gamma_coincidence
                    new_data_list.append(new_data_point)
            else:
                # This data point has a different radionuclide
                pass
    data_list = new_data_list

    return data_list


# NOTE OK ABOVE!

def analyze_files(data_list, plotcard, metadata, data_path):
    print("Starting analysis of data files!")

    # Need to iterate through the relevant files
    for i, (key, value) in enumerate(metadata.items()):

        # NOTE just for testing
        if i > 4:
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

        # Based on the file metadata, this should be the data point
        data_point_criteria = {
            "detector_type": value["properties"]["detector_type"],
            "detector_diameter": value["properties"]["detector_diameter"],
            "detector_length": value["properties"]["detector_length"],
            "detector_source_distance": value["properties"]["detector_source_distance"],
            "source_type": value["properties"]["source_type"],
        }

        if run_type == "radionuclides":
            # Based on the file, this should be the ZA of the data_point
            ZA_criteria = [value["properties"]["Z"], value["properties"]["A"]]
            data_point_criteria["ZA"] = ZA_criteria

            for j, ZA in enumerate(plotcard["radionuclides"]["ZAs"]):
                if ZA == ZA_criteria:
                    # Check which energies are relevant for this radionuclide
                    gammas_singles = plotcard["radionuclides"]["gammas_singles"][j]
                    gammas_coincidences = plotcard["radionuclides"]["gammas_coincidences"][j]

                    # Analyze single gamma rays
                    for E_gamma in gammas_singles:
                        if E_gamma == "combined":
                            # Handle this later
                            counts = "n/a"
                            effint, effint_unc = 0.0, 0.0
                        else:
                            ROI_polynial_coeffs = plotcard["analysis"]["ROI_width_polynomial"]
                            ROI_width = np.polyval(ROI_polynial_coeffs, E_gamma)

                            counts_a = ROI_analysis_1D(events_single_a, E_gamma, ROI_width)
                            counts_b = ROI_analysis_1D(events_single_b, E_gamma, ROI_width)
                            print("\t\t\tGamma ray: " + str(E_gamma) + " keV, counts a: " + str(counts_a) + ", counts b: " + str(counts_b))

                            counts = counts_a + counts_b
                            events = value["properties"]["events"]

                            # Caluclate the effint
                            effint, effint_unc = caluclate_effint(counts, events)

                        this_data_point_criteria = data_point_criteria.copy()
                        this_data_point_criteria["analysis_type"] = "singles"
                        this_data_point_criteria["Egamma"] = E_gamma

                        # Store the data in the data_list
                        for data_point in data_list:
                            if all(data_point[k] == v for k, v in this_data_point_criteria.items()):
                                data_point["counts_radionuclide"] = counts
                                data_point["effint"] = [effint, effint_unc]

                    # Analyze coincidence gamma rays also
                    for Egamma1_Egamma2 in gammas_coincidences:
                        if Egamma1_Egamma2 == "combined":
                            # Handle this later
                            counts = "n/a"
                            effint, effint_unc = 0.0, 0.0
                        else:
                            E_gamma1, E_gamma2 = Egamma1_Egamma2

                            ROI_polynial_coeffs = plotcard["analysis"]["ROI_width_polynomial"]
                            ROI_width1 = np.polyval(ROI_polynial_coeffs, E_gamma1)
                            ROI_width2 = np.polyval(ROI_polynial_coeffs, E_gamma2)

                            counts_a1b2 = ROI_analysis_2D(events_coincidence_a, events_coincidence_b, E_gamma1, E_gamma2, ROI_width1, ROI_width2)
                            counts_b1a2 = ROI_analysis_2D(events_coincidence_b, events_coincidence_a, E_gamma1, E_gamma2, ROI_width1, ROI_width2)
                            print("\t\t\tGamma ray 1: " + str(E_gamma1) + " keV, gamma ray 2: " + str(E_gamma2) + " keV , counts a1b2: " + str(counts_a1b2) + ", counts b2a1: " + str(counts_b1a2))

                            counts = counts_a1b2 + counts_b1a2
                            events = value["properties"]["events"]

                            # Caluclate the effint
                            effint, effint_unc = caluclate_effint(counts, events)

                        this_data_point_criteria = data_point_criteria.copy()
                        this_data_point_criteria["analysis_type"] = "coincidences"
                        this_data_point_criteria["Egamma1_Egamma2"] = Egamma1_Egamma2

                        # Store the data in the data_list
                        for data_point in data_list:
                            if all(data_point[k] == v for k, v in this_data_point_criteria.items()):
                                data_point["counts_radionuclide"] = counts
                                data_point["effint"] = [effint, effint_unc]

        elif run_type == "background":

            # for j, ZA in enumerate(plotcard["radionuclides"]["ZAs"]):
            #         print("\t\tAnalyzing radionuclide " + str(j+1) + " out of " + str(len(plotcard["radionuclides"]["ZAs"])) + f": Z={ZA[0]}, A={A[1]}")

            # iterate through all the radionuclides and gamma rays here (all of them!)

            pass
        elif run_type == "filter":
            # Do this later
            pass
        else:
            pass

                


                

        #         counts_coincidences.append([counts_a1b2, counts_b1a2])
            
        #     radionuclide_counts_singles.append(counts_singles)
        #     radionuclide_counts_coincidences.append(counts_coincidences)
        
        # value["counts_singles"] = radionuclide_counts_singles
        # value["counts_coincidences"] = radionuclide_counts_coincidences
        # plotcard["metadata"][key] = value

    return data_list


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


    # TODO save as list of dictionaries instead!!!
    # each dictionary like
    #{"detector_type": X, "detector_diameter": X, "detector_length": X, "detector_source_distance": X, "source_type": X, "ZA": X, "analysis_type": singles/coincidences, "gammas_XXX": X, "background_file": X,       NOW THE DATA PART: "" "background counts"}

    # pprint.pp(plotcard)

    return plotcard


def caluclate_effint(total_counts, events):
    # Calculate efficiency and intensity
    effint = total_counts / events

    # Use binomial instead of poisson approximation
    effint_unc = np.sqrt(total_counts * (1 - total_counts/events)) / events

    return float(effint), float(effint_unc)


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
            

def plot_mda(plotcard):
    metadata = plotcard[""]

    # Now the MDA can be caluclated in the next step
    return


def calculate_mda(LD, effint, tM, t12=0):
    # Calculate the minimum detectable activity
    mda = LD / (effint * tM)

    # Optional correction for decay during measurement
    if t12 != 0:
        lambdaa = np.log(2)/t12 # extra a to respect Python reserved word
        decay_correction = (lambdaa * tM) / (1 - np.exp(-lambdaa * tM))
        mda *= decay_correction

    return mda


def calculate_combined_mda():

    # Check if the combined MDA should be calculated
    # if "combined" in plotcard["radionuclides"]["gammas_singles"][j]:
    #     print("do combined also")
    # This should be for every radionucldie... need to iterate through ZAs
    return


def make_plot():
    return


# step 3: mda calculation

# step 4: plot and save








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

for data_point in data_list:
    print(data_point["effint"])


# plotcard = calculate_quantities(plotcard)

# plot_mda(plotcard)

# with open("temp.json", "w") as f:
#     json.dump(plotcard, f, indent=4)

# print(plotcard)

# print(metadata)


