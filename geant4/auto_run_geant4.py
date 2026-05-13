import argparse
import json
import time




def run_geometry(geometry, Z=0, A=0):
    detector_type = geometry["detector_type"]
    detector_diameter = geometry["detector_diameter"]
    detector_length = geometry["detector_length"]
    detector_source_distance = geometry["detector_source_distance"]

    if len(detector_length) == 1 and type(detector_length[0]) == str:
        detector_length = geometry[detector_length[0]]

    # Combine all geometry to iterate over
    settings = [(det_type, det_diam, det_leng, det_sdis) 
                for det_type in detector_type
                for det_diam in detector_diameter
                for det_leng in detector_length
                for det_sdis in detector_source_distance]

    for i_s, (det_type, det_diam, det_leng, det_sdis) in enumerate(settings):
        print(f"\t\tRunning: detector_type={det_type}, detector_diameter={det_diam}, detector_length={det_leng}, detector_source_distance={det_sdis}")
        # print(det_type, det_diam, det_leng, det_sdis)

    # # Iterate over radionuclides
    # for i_s, (ZA, diameter, length, distance) in enumerate(settings):
    #     # Optional: skip some iterations
    #     # if i_s <= 2:
    #     #     continue
    #     pass


def run_radionuclides(radionuclides, geometry):
    if radionuclides["active"] == False:
        return
    
    ZAs = radionuclides["ZAs"]
    for Z, A in ZAs:
        print(f"\tRunning: Z={Z}, A={A}")

        run_geometry(geometry, Z=Z, A=A)




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

