import numpy as np
import uproot
import argparse
# from numba import jit
import json
import re
# from uncertainties import ufloat


def start_analysis(ZAs):
    print("Starting analysis: ")
    for Z, A in ZAs:
        print("\tZ=" + str(Z) + ", A=" + str(A))

        # Go through the metadata and look for the interesting simulations
        data_keys = []
        data_filenames = []
        for key, value in metadata.items():
            if value["type"] != "radionuclides":
                continue
            if value["properties"]["Z"] != Z or value["properties"]["A"] != A:
                continue
            data_filenames.append(value["filename"])
            data_keys.append(key)
        
        print("\tFound " + str(len(data_filenames)) + " files")

        # Now do the analysis
        for i, data_filename in enumerate(data_filenames):
            print("\t\tAnalyzing file " + str(i+1) + " out of " + str(len(data_filenames)) + ": " + str(data_filename))

            analyze_single_and_coincidence(Z, A, data_filename)


def analyze_single_and_coincidence(Z, A, data_filename):
    data_filepath = args.data + data_filename

    data_file = uproot.open(data_filepath)
    tree_detector_a = data_file["E_a_list"]
    tree_detector_b = data_file["E_b_list"]
    tree_detector_ab = data_file["E_ab_list"]

    # Single gamma-rays
    events_single_a = tree_detector_a["energy_a"].array(library="np") * 1e3
    events_single_b = tree_detector_b["energy_b"].array(library="np") * 1e3

    # Coincidence gamma-rays
    events_coincidence_a = tree_detector_ab["energy_a"].array(library="np") * 1e3
    events_coincidence_b = tree_detector_ab["energy_b"].array(library="np") * 1e3

    # Gamma ray energy infromation
    gammas_single = gammas["Z" + str(Z) + "_A" + str(A)]["singles"]["E_gamma"]
    gammas_coincidence_1 = gammas["Z" + str(Z) + "_A" + str(A)]["coincidences"]["E_gamma1"]
    gammas_coincidence_2 = gammas["Z" + str(Z) + "_A" + str(A)]["coincidences"]["E_gamma2"]

    # Single detector counts
    counts_single_a = []
    counts_single_b = []
    for E_gamma in gammas_single:
        counts_a = count_1d(events_single_a, E_gamma, ROI_standard)
        counts_single_a.append(int(counts_a))
        counts_b = count_1d(events_single_b, E_gamma, ROI_standard)
        counts_single_b.append(int(counts_b))

    # Coincidence counts
    counts_coincidence_a1b2 = []
    counts_coincidence_a2b1 = []
    for E_gamma1, E_gamma2 in zip(gammas_coincidence_1, gammas_coincidence_2):
        counts_a1b2 = count_2d(events_coincidence_a, events_coincidence_b, E_gamma1, E_gamma2, ROI_standard)
        counts_coincidence_a1b2.append(int(counts_a1b2))
        counts_a2b1 = count_2d(events_coincidence_a, events_coincidence_b, E_gamma2, E_gamma1, ROI_standard)
        counts_coincidence_a2b1.append(int(counts_a2b1))

    # print(gammas_single)
    # print(counts_single_a)
    # print(counts_single_b)
    # print(gammas_coincidence_1)
    # print(gammas_coincidence_2)
    # print(counts_coincidence_a1b2)
    # print(counts_coincidence_a2b1)
    return


# @jit(nopython=True)
def count_1d(events, Eg, ROI):
    # Energy distance is half of the ROI size
    dE = ROI/2

    # Count peak counts
    cond = np.logical_and(events > Eg-dE, events < Eg+dE)
    counts = cond.sum()
    return counts


def count_2d(events_a, events_b, Eg1, Eg2, ROI):
    # Energy distance is half of the ROI size
    dE = ROI/2

    # Count peak counts
    cond_1 = np.logical_and(events_a > Eg1-dE, events_a < Eg1+dE)
    cond_2 = np.logical_and(events_b > Eg2-dE, events_b < Eg2+dE)
    cond = np.logical_and(cond_1, cond_2)
    counts = cond.sum()
    return counts


# Parser for adding arguments
parser = argparse.ArgumentParser(prog="DetermineEffInt",
                                 description="Analyze .root files to determine the product of efficiency and intensity",
                                 epilog="Elias Arnqvist, 2026, Uppsala University",
                                 add_help=True)
parser.add_argument("-Z", type=int, required=False, default=0, help="Atomic number (Z) of radionuclide")
parser.add_argument("-A", type=int, required=False, default=0, help="Mass number (A) of radionuclide")
parser.add_argument("-d", "--data", type=str, required=False, default="../../geant4/output_pelle/", help="Path to directory with root files")
parser.add_argument("-m", "--metadata", type=str, required=False, default="../../geant4/output_pelle/metadata.json", help="Path to metadata file")
parser.add_argument("-g", "--gammas", type=str, required=False, default="gammas.json", help="Path to gamma ray information")
parser.add_argument("-s", "--save", type=str, required=False, default="/efficiency_intensity.json", help="Path to save output data")
args = parser.parse_args()

ROI_standard = 5

# Open gamma ray data file
gammas_filepath = args.gammas
print("Opening gamma file: " + gammas_filepath)
with open(gammas_filepath, "r") as f:
    gammas = json.load(f)

# Open metadata file
metadata_filepath = args.metadata
print("Opening metadata file: " + metadata_filepath)
with open(metadata_filepath, "r") as f:
    metadata = json.load(f)

# Dictionary to save results in
output = {}

# Look for radionuclides in the gammas and metadata
# Make a list of Z, A present in gammas and metadata
ZAs_gammas = [(int(m.group(1)), int(m.group(2))) for key in gammas.keys() if (m := re.match(r"Z(\d+)_A(\d+)", key))]
# Automatically no duplicates here from json load
ZAs_metadata = [(value["properties"]["Z"], value["properties"]["A"]) for value in metadata.values() if value["type"] == "radionuclides"]
ZAs_metadata = list(set(ZAs_metadata))
# If Z and A are specified use them, otherwise use all radionuclides
Z = args.Z
A = args.A
if Z != 0 and A != 0:
    # Z and A are specified
    print("Using specified ZA: Z=" + str(Z) + ", A=" + str(A))
    if (Z, A) not in ZAs_gammas:
        raise Exception("Used Z and A not found in gamma ray data")
    if (Z, A) not in ZAs_metadata:
        raise Exception("Used Z and A not found in metadata")

    start_analysis([(Z, A)])
else:
    # Z and A not specified, so use the ones in the gammas
    print("Using ZA in gamma file: " + str(ZAs_gammas))
    for ZA in ZAs_gammas:
        if ZA not in ZAs_metadata:
            print("Could not find " + str(ZA) + " in metadata!")
            raise Exception("Used Z and A not found in metadata")

