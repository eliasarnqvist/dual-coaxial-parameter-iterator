import argparse
import json
import time
import subprocess
import uuid
import os
import numpy as np


def run_geometry(geometry, run_dict, run_type):
    detector_type = geometry["detector_type"]
    detector_diameter = geometry["detector_diameter"]
    detector_length = geometry["detector_length"]
    detector_source_distance = geometry["detector_source_distance"]
    detector_source_type = geometry["source_type"]

    # If this is a background run, prepare the SURE source spectrum
    if run_type == "background":
        # Background source spectrum for SURE model
        background_filename = run_dict["source_spectrum"]
        flux_data = np.genfromtxt(background_filename, delimiter=",")
        E = flux_data[:, 0]
        # Flux is in s^(-1)cm^(-2) here (not mm^(-2))
        flux = flux_data[:, 1]
        background_total_flux = flux.sum()
        # Make cumulative distribution function of the flux
        flux_cdf = np.cumsum(flux) / background_total_flux
        # Save with all 9 decimals
        np.savetxt("resources/flux_cdf.dat", np.column_stack((E, flux_cdf)), fmt="%.9f\t%.9f")

    # Combine all geometry to iterate over
    settings = [(det_type, det_diam, det_leng, det_sdis, det_styp) 
                for det_type in detector_type
                for det_diam in detector_diameter
                for det_leng in detector_length
                for det_sdis in detector_source_distance
                for det_styp in detector_source_type]

    for i_s, (det_type, det_diam, det_leng, det_sdis, det_styp) in enumerate(settings):
        # Special case to put the length and diameter to the same value
        if det_leng == "detector_diameter":
            det_leng = det_diam
        else:
            print("Inclorrect detector diameter specified!")
            break

        # Optional: skip some iterations
        # if i_s <= 2:
        #     continue

        print(f"\t\tRunning: detector_type={det_type}, detector_diameter={det_diam}, detector_length={det_leng}, detector_source_distance={det_sdis}, detector_source_type={det_styp}")

        # Make macro file first
        macro_content = ""

        # Threads and filename
        threads = run_dict["threads"]
        macro_content += "/run/numberOfThreads " + str(threads) + "\n"
        file_name = run_dict["output"] + "threadoutput_" + str(i_s) + "_.root"
        macro_content += "/E_file_settings/fileName " + file_name + "\n"

        # Dtector geometry settings
        macro_content += "/E_detector/detectorDiameter " + str(det_diam) + "\n"
        macro_content += "/E_detector/detectorLength " + str(det_leng) + "\n"
        macro_content += "/E_detector/sourceDistance " + str(det_sdis) + "\n"
        macro_content += "/E_detector/detectorType " + str(det_type) + "\n"
        macro_content += "/E_detector/sourceType " + str(det_styp) + "\n"

        # General settings
        macro_content += "/run/reinitializeGeometry" + "\n"
        macro_content += "/run/initialize" + "\n"
        macro_content += "/process/had/rdm/thresholdForVeryLongDecayTime 1.0e+60 year" + "\n"

        # If this is a radionuclide run
        if run_type == "radionuclides" or run_type == "filter":
            # Specify the ion
            macro_content += "/gun/particle ion" + "\n"
            Z, A = run_dict["ZA"]
            macro_content += "/gun/ion " + str(Z) + " " + str(A) + " 0 0" + "\n"
            macro_content += "/process/had/rdm/nucleusLimits "+str(A)+" "+str(A)+" "+str(Z)+" "+str(Z)+"\n"

            # Specify the source type (either point or filter - no SURE model here)
            macro_content += "/E_source/sourceType " + str(det_styp) + "\n"
        elif run_type == "background":
            # The background model uses the SURE model
            # Calculate the SURE model radius that encompasses both detectors
            SURE_radius = calculate_SURE_radius(det_diam, det_leng, det_styp, det_type)

            # The SURE model must be used (source type 2) with a specified SURE radius
            macro_content += "/E_source/sourceType " + str(2) + "\n"
            macro_content += "/E_source/sourceRadiusSURE " + str(SURE_radius) + "\n"

        # General run settings
        events = run_dict["events"]
        macro_content += "/run/printProgress " + str(int(events/10)) + "\n"
        macro_content += "/run/beamOn " + str(int(events))

        # Save the macro file now
        print("\t\t\tWriting macro file...")
        build_folder = "build/"
        macro_name = "autorun.mac"
        with open(build_folder + macro_name, "w") as file:
            file.write(macro_content)

        # Make the metadata file if it does not exist already
        metadata_folder = run_dict["metadata"]
        if not os.path.exists(metadata_folder + "metadata.json"):
            with open(metadata_folder + "metadata.json", "w") as f:
                f.write("{}")
        
        # Open the metadata file
        with open(metadata_folder + "metadata.json") as f:
            metadata = json.load(f)
        
        # Check if this run already exists in the metadata
        already_simulated = any(
            present_entry["type"] == run_type and
            present_entry["properties"]["detector_type"] == det_type and
            present_entry["properties"]["detector_diameter"] == det_diam and
            present_entry["properties"]["detector_length"] == det_leng and
            present_entry["properties"]["source_distance"] == det_sdis and
            present_entry["properties"]["source_type"] == det_styp and
            ((present_entry["properties"]["Z"] == Z) if (run_type == "radionuclides" or run_type == "filter") else True) and
            ((present_entry["properties"]["A"] == A) if (run_type == "radionuclides" or run_type == "filter") else True) and
            ((present_entry["properties"]["background_file"] == background_filename) if (run_type == "background") else True)
            for present_entry in metadata.values()
        )

        if not already_simulated and not run_dict["test_run"]:
            sim_start_time = time.time()

            # Start the Geant4 simulation
            print("\t\t\tRunning Geant4...")
            process_geant4 = [build_folder + "sim", build_folder + macro_name]
            # Top option is without verbocity, bottom is with verbocity
            result = subprocess.run(process_geant4, stdout=subprocess.DEVNULL)
            # result = subprocess.run(process_geant4) # DO NOT PUT SHELL=True

            sim_stop_time = time.time()
            # The time it took to only run geant4
            simulated_minutes = (sim_stop_time - sim_start_time) / 60

            # Combine the ROOT files and give it a random and unique uuid4 name
            print("\t\t\tCombining ROOT files...")
            run_id = str(uuid.uuid4())
            output_folder = run_dict["output"]
            output_file = output_folder + run_id + ".root"
            # Combine thread output files with the built-in ROOT hadd
            process_root = "hadd -f " + output_file + " " + output_folder + "threadoutput_" + str(i_s) + "*.root"
            result = subprocess.run(process_root, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            # result = subprocess.run(process_root, shell=True)

            # Now add the run metadata to the metadata file
            print("\t\t\tAdding metadata...")

            # Add the run information to the metadata
            # Properties will depend on what type of run this is
            properties = {
                # Properties related to the detector and source type
                "detector_type":int(det_type),
                "detector_diameter":det_diam,
                "detector_length":det_leng,
                "source_distance":det_sdis,
                "source_type":int(det_styp),
                # Properties related to the simulation
                "events":int(events),
                "threads":int(threads),
                "time_minutes":simulated_minutes,
                "throughput":events/(simulated_minutes*60*threads),
                }
            if run_type == "radionuclides" or run_type == "filter":
                # Save the Z and A of the simulated radionuclide (same for filter)
                Z, A = run_dict["ZA"]
                properties["Z"] = int(Z)
                properties["A"] = int(A)
            elif run_type == "background":
                # Save the properties of the SURE background model
                properties["background_file"] = background_filename
                properties["SURE_radius"] = SURE_radius

                # Equivalent real time of the simulation, also known as pseudo time
                # Note flux conversion from cm^-2 to mm^-2 to be compatible with sure radius in mm
                # Or an equivalent would be to convert radius from mm to cm here
                pseudo_time = events / (SURE_radius**2 * np.pi * (background_total_flux / 100))
                properties["SURE_pseudo_time"] = pseudo_time
            
            # Now we can put together the metadata
            metadata[run_id] = {
                # General metadata first
                "filename":(run_id + ".root"),
                "file_size":os.path.getsize(output_file),
                "type":run_type,
                # Properties from above
                "properties":properties
                }
            
            # Write the updated metadata
            with open(metadata_folder + "metadata.json", "w") as f:
                json.dump(metadata, f, indent=4)

            print("\t\t\tDeleting temporary ROOT files...")
            process_delete = "rm " + output_folder + "threadoutput_" + str(i_s) + "*.root"
            # result = subprocess.run(process_delete, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            result = subprocess.run(process_delete, shell=True)

            print(f"\t\t\tTime spent for previous run: {simulated_minutes:.2f} minutes")
        elif not already_simulated and run_dict["test_run"]:
            print(f"\t\t\tWould have simulated this! But test_run = true")
        else:
            print(f"\t\t\tFound this run in metadata already, skipping!")


def calculate_SURE_radius(diameter, length, source_type, detector_type):
    # Correct for filter source thickness if needed
    if source_type == 1:
        # Filter source
        source_thickness = 13.8/2
    elif source_type == 0:
        # Point source
        source_thickness = 0
    
    # Determine radius differently for coaxial and planar detectors
    if detector_type == 0 or detector_type == 1:
        # Coaxial p- or n-type
        R = np.sqrt((length+1.5+4+5+5+5+5+source_thickness)**2 + (diameter/2+2+4+1.5)**2)
    elif detector_type == 2:
        # Planar
        R = np.sqrt((length+1.5+4+5+2+5+5+source_thickness)**2 + (diameter/2+2+4+1.5)**2)
    
    return R


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
    if filter["active"] == False:
        return
    else:
        ZAs = filter["ZAs"]
        for i_ZA, (Z, A) in enumerate(ZAs):
            print(f"\tRunning: Z={Z}, A={A} ({i_ZA+1} out of {len(ZAs)})")

            filter["ZA"] = [Z, A]
            run_geometry(geometry, run_dict=filter, run_type="filter")


def run(runcard):
    # Keep track of the time the run takes
    start_time = time.time()

    if not "geometry" in runcard.keys():
        print("Please specify the geometry in the runcard")
        return

    if "radionuclides" in runcard.keys():
        print("Running radionucldies")
        run_radionuclides(runcard["radionuclides"], runcard["geometry"])
    if "background" in runcard.keys():
        print("Running background")
        run_background(runcard["background"], runcard["geometry"])
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

