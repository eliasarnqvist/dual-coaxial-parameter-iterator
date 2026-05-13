import argparse
import json
import time

def iterate_geometry():
    
    settings = [(ZA, dia, len, dis) 
                for ZA in ZAs
                for dia in detector_diameters
                for len in detector_lengths
                for dis in source_distances]

    # Iterate over radionuclides
    for i_s, (ZA, diameter, length, distance) in enumerate(settings):
        # Optional: skip some iterations
        # if i_s <= 2:
        #     continue
        pass


def run_radionuclides(radionuclides, geometry):
    if radionuclides["active"] == False:
        return
    
    ZAs = radionuclides["ZAs"]


def run_background(background, geometry):
    return


def run_filter(filter, geometry):
    return


def run(runcard):
    # Keep track of the time the run takes
    start_time = time.time()
    time_interval = start_time

    if not "geometry" in runcard.keys():
        print("Please specify the geometry in the runcard")
        return

    if "radionuclides" in runcard.keys():
        print("Running radionucldies")
        run_radionuclides(runcard["radionuclides"], runcard["geometry"])
    if "background" in runcard.keys():
        print("Running background")
        run_background(runcard["background", runcard["geometry"]])
    if "filter" in runcard.keys():
        print("Running filter background")
        run_filter(runcard["filter"], runcard["geometry"])

    

    print("Finished!")

    # Print the total time this runcard took
    end_time = time.time()
    elapsed_minutes = (end_time - start_time) / 60
    print(f"Total time spent: {elapsed_minutes:.2f} minutes")



# Parser for adding arguments
parser = argparse.ArgumentParser(prog="auto_run_geant4",
                                 description="Run geant4 and save data",
                                 epilog="Elias Arnqvist, 2026, Uppsala University",
                                 add_help=True)
parser.add_argument("-r", "--runcard", type=str, required=True, default="run.json", help="File specifying what to run")
parser.add_argument("-s" "--save_folder", type=str, required=False, default="output", help="Output folder")
args = parser.parse_args()

runcard_filepath = args.runcard
with open(runcard_filepath, "r") as f:
    runcard = json.load(f)

print(runcard)

run(runcard)

