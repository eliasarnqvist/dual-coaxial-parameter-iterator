import argparse
import json
import time
import subprocess
import uuid
import os


def run_geometry(geometry, run_dict, run_type):
    detector_type = geometry["detector_type"]
    detector_diameter = geometry["detector_diameter"]
    detector_length = geometry["detector_length"]
    detector_source_distance = geometry["detector_source_distance"]
    detector_source_type = geometry["source_type"]

    # Special case to put the length and diameter to the same value
    if len(detector_length) == 1 and type(detector_length[0]) == str:
        detector_length = geometry[detector_length[0]]

    # Combine all geometry to iterate over
    settings = [(det_type, det_diam, det_leng, det_sdis, det_styp) 
                for det_type in detector_type
                for det_diam in detector_diameter
                for det_leng in detector_length
                for det_sdis in detector_source_distance
                for det_styp in detector_source_type]

    for i_s, (det_type, det_diam, det_leng, det_sdis, det_styp) in enumerate(settings):
        # Optional: skip some iterations
        # if i_s <= 2:
        #     continue

        print(f"\t\tRunning: detector_type={det_type}, detector_diameter={det_diam}, detector_length={det_leng}, detector_source_distance={det_sdis}, detector_source_type={det_styp}")
            
        # Make macro file first
        macro_content = ""

        # Threads and filename
        threads = run_dict["threads"]
        macro_content += "/run/numberOfThreads " + str(threads) + "\n"
        file_name = run_dict["output"] + "threadoutput_" + str(i_s)
        macro_content += "/E_file_settings/fileName " + file_name + "\n"

        # Dtector geometry settings
        macro_content += "/E_detector/detectorDiameter " + str(det_diam) + "\n"
        macro_content += "/E_detector/detectorLength " + str(det_leng) + "\n"
        macro_content += "/E_detector/sourceDistance " + str(det_sdis) + "\n"
        macro_content += "/E_detector/detectorType " + str(det_type) + "\n" # TODO MUST IMPLEMENT THIS in G4
        macro_content += "/E_detector/sourceType " + str(det_styp) + "\n" # TODO MUST IMPLEMENT THIS in G4

        # General settings
        macro_content += "/run/reinitializeGeometry" + "\n"
        macro_content += "/run/initialize" + "\n"
        macro_content += "/process/had/rdm/thresholdForVeryLongDecayTime 1.0e+60 year" + "\n"

        # If this is a radionuclide run
        if run_type == "radionuclides":
            # Specify the ion
            macro_content += "/gun/particle ion" + "\n"
            Z, A = run_dict["ZA"]
            macro_content += "/gun/ion " + str(Z) + " " + str(A) + " 0 0" + "\n"
            macro_content += "/process/had/rdm/nucleusLimits "+str(A)+" "+str(A)+" "+str(Z)+" "+str(Z)+"\n"

            # Specify the source type (either point or filter - no SURE model here)
            macro_content += "/E_source/sourceType " + str(det_styp) + "\n"
        elif run_type == "background":

            #TODO Implement sure radius calculation depending on detector type and size

            # "source_type":source_type,
            # "source_SURE_radius":source_SURE_radius,
            # "SURE_background_total_flux":background_total_flux/100, # convert to mm^-2 s^-1
            # "SURE_pseudo_time":pseudo_time,

            # The sure model must be used with a specified SURE radius
            macro_content += "/E_source/sourceType " + str(det_styp) + "\n"
            macro_content += "/E_source/sourceRadiusSURE " + 1 + "\n"

        # General run settings
        events = run_dict["events"]
        macro_content += "/run/printProgress " + str(int(events/10)) + "\n"
        macro_content += "/run/beamOn " + str(int(events))

        # Save the macro file now
        print("\tWriting macro file...")
        build_folder = "build/"
        macro_name = "autorun.mac"
        with open(build_folder + macro_name, "w") as file:
            file.write(macro_content)

        sim_start_time = time.time()

        # Start the Geant4 simulation
        print("\tRunning Geant4...")
        process_geant4 = [build_folder + "sim", build_folder + macro_name]
        # Top option is without verbocity, bottom is with verbocity
        result = subprocess.run(process_geant4, stdout=subprocess.DEVNULL)
        # result = subprocess.run(process_geant4) # DO NOT PUT SHELL=True

        sim_stop_time = time.time()
        # The time it took to only run geant4
        simulated_minutes = (sim_stop_time - sim_start_time) / 60

        # Combine the ROOT files and give it a random uuid4 name
        print("\tCombining ROOT files...")
        run_id = str(uuid.uuid4())
        output_folder = run_dict["output"]
        output_file = output_folder + run_id + ".root"
        # Combine with ROOT hadd
        process_root = "hadd -f " + output_file + " " + output_folder + "threadoutput_" + str(i_s) + "*.root"
        result = subprocess.run(process_root, shell=True, stdout=subprocess.DEVNULL)
        # result = subprocess.run(process_root, shell=True)

        # Now add the run metadata to the metadata file
        print("\tAdding metadata...")
        # Make the file if it does not exist already
        if not os.path.exists(output_folder + "metadata.json"):
            with open(output_folder + "metadata.json", "w") as f:
                f.write("{}")
        # Open the metadata file
        with open(output_folder + "metadata.json") as f:
            metadata = json.load(f)
        # Add the run information to the metadata
        # Properties will depend on what type of run this is
        properties = {
            # Properties related to the detector and source type
            "detector_type":det_type,
            "detector_diameter":det_diam,
            "detector_length":det_leng,
            "source_distance":det_sdis,
            "source_type":det_styp,
            # Properties related to the simulation
            "events":events,
            "threads":threads,
            "time_minutes":simulated_minutes,
            "throughput":events/(simulated_minutes*60*threads),}
        if run_type == "radionuclides":
            Z, A = run_dict["ZA"]
            properties["Z"] = Z
            properties["A"] = A
        elif run_type == "background":
            properties["background_file"] = run_dict["source_spectrum"]

            # TODO 
            # "source_SURE_radius":source_SURE_radius,
            # "SURE_background_total_flux":background_total_flux/100, # convert to mm^-2 s^-1
            # "SURE_pseudo_time":pseudo_time,

        elif run_type == "filter":
            Z, A = run_dict["ZA"]
            properties["Z"] = Z
            properties["A"] = A

        metadata[run_id] = {
            # General metadata first
            "filename":(run_id + ".root"),
            "file_size":os.path.getsize(output_file),
            "type":run_type,
            # Properties from above
            "properties":properties
            }
        # Write the updated metadata
        with open(output_folder + "metadata.json", "w") as f:
            json.dump(metadata, f, indent=4)

        print("\tDeleting temporary ROOT files...")
        process_delete = "rm " + output_folder + "threadoutput_" + str(i_s) + "*.root"
        # result = subprocess.run(process_delete, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        result = subprocess.run(process_delete, shell=True)

        partial_time = time.time()
        elapsed_minutes = (partial_time - time_interval) / 60
        print(f"\tTime spent for previous run: {elapsed_minutes:.2f} minutes")
        time_interval = partial_time


def run_radionuclides(radionuclides, geometry):
    if radionuclides["active"] == False:
        return
    else:
        ZAs = radionuclides["ZAs"]
        for i_ZA, (Z, A) in enumerate(ZAs):
            print(f"\tRunning: Z={Z}, A={A} ({i_ZA+1} out of {len(ZAs)})")

            radionuclides["ZA"] = [Z, A]
            run_geometry(geometry, run_dict=radionuclides, run_type="radionuclides")


def run_background(background, geometry):
    if background["active"] == False:
        return
    else:
        source_spectra = background["source_spectra"]
        for i_source_spectrum, source_spectrum in enumerate(source_spectra):
            print(f"\tRunning: source_spectrum={source_spectrum} ({i_source_spectrum+1} out of {len(source_spectra)})")

            background["source_spectrum"] = source_spectrum
            run_geometry(geometry, run_dict=background, run_type="background")


def run_filter(filter, geometry):
    return
# TODO fix this function


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

