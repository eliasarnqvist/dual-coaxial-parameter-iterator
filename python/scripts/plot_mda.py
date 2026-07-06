import argparse
import json
import time
import subprocess
import uuid
import os
import numpy as np


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

def analyze_file(plotcard, metadata, data_filepath):

    for i, (key, value) in enumerate(metadata.items()):
        data_filename = value["filename"]
        data_filesize = value["file_size"] / (1024**2)
        print("\tAnalyzing file " + str(i+1) + " out of " + str(len(metadata)) + ": " + str(data_filename) + f" ({data_filesize:.2f} mb)")


    # result should be LD, effint, mda, x value, and line label


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
data_filepath = args.data

metadata = cut_metadata(plotcard, metadata)

analyze_file(plotcard, metadata, data_filepath)


# print(plotcard)

# print(metadata)


