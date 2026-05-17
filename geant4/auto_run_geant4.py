import argparse
import json
import time


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
        threads = run_dict["threads"]
        macro_content += "/run/numberOfThreads " + str(threads) + "\n"

        file_name = run_dict["output"] + "threadoutput_" + run_type + "_" + str(i_s)
        macro_content += "/E_file_settings/fileName " + file_name + "\n"

        macro_content += "/E_detector/detectorDiameter " + str(det_diam) + "\n"
        macro_content += "/E_detector/detectorLength " + str(det_leng) + "\n"
        macro_content += "/E_detector/sourceDistance " + str(det_sdis) + "\n"
        macro_content += "/E_detector/detectorType " + str(det_type) + "\n" # MUST IMPLEMENT THIS in G4
        macro_content += "/E_detector/sourceType " + str(det_styp) + "\n" # MUST IMPLEMENT THIS in G4

        macro_content += "/run/reinitializeGeometry" + "\n"
        macro_content += "/run/initialize" + "\n"
        macro_content += "/process/had/rdm/thresholdForVeryLongDecayTime 1.0e+60 year" + "\n"

        if run_type == "radionuclides":
            macro_content += "/gun/particle ion" + "\n"
            Z, A = run_dict["ZA"]
            macro_content += "/gun/ion " + str(Z) + " " + str(A) + " 0 0" + "\n"
            macro_content += "/process/had/rdm/nucleusLimits "+str(A)+" "+str(A)+" "+str(Z)+" "+str(Z)+"\n"

        macro_content += "/E_source/sourceType " + str(source_type) + "\n"
        # macro_content += "/E_source/sourceRadiusSURE " + XXX + "\n"

        macro_content += "/run/printProgress " + str(int(events_per_radionuclide/10)) + "\n"
        macro_content += "/run/beamOn " + str(int(events_per_radionuclide))

        print("\tWriting macro file...")

        with open(build_folder + macro_name, "w") as file:
            file.write(macro_content)

        sim_start_time = time.time()

        print("\tRunning Geant4...")
        process_geant4 = [build_folder + "sim", build_folder + macro_name]
        result = subprocess.run(process_geant4, stdout=subprocess.DEVNULL)
        # result = subprocess.run(process_geant4) # DO NOT PUT SHELL=True

        sim_stop_time = time.time()
        # the time it took to only run geant4
        simulated_minutes = (sim_stop_time - sim_start_time) / 60

        print("\tCombining ROOT files...")
        run_id = str(uuid.uuid4())
        # output_file = output_folder + "nucOutput_" + str(i_s) + settings_name
        output_file = output_folder + run_id + ".root"
        process_root = "hadd -f " + output_file + " " + output_folder + "threadoutput_" + str(i_s) + "*.root"
        result = subprocess.run(process_root, shell=True, stdout=subprocess.DEVNULL)
        # result = subprocess.run(process_root, shell=True)

        print("\tAdding metadata...")
        if not os.path.exists(output_folder + "metadata.json"):
            with open(output_folder + "metadata.json", "w") as f:
                f.write("{}")
        with open(output_folder + "metadata.json") as f:
            metadata = json.load(f)
        metadata[run_id] = {
            "filename":(run_id + ".root"),
            "file_size":os.path.getsize(output_file),
            "type":"radionuclides",
            "properties":{
                "model":"coaxial_v2",
                "Z":Z,
                "A":A,
                "detector_diameter":diameter,
                "detector_length":length,
                "source_distance":distance,
                "select_ntype_instead_of_ptype":select_n_type_instead_of_p_type,
                "select_filter_source":select_filter_source,
                "source_type":source_type,
                "source_SURE_radius":source_SURE_radius,
                "source_type":"FOI_filter_v1",
                "events":events_per_radionuclide,
                "threads":number_of_threads,
                "time":simulated_minutes,
                "throughput":events_per_radionuclide/(simulated_minutes*60*number_of_threads),
                },
            }
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
    
    ZAs = radionuclides["ZAs"]
    for i_ZA, (Z, A) in enumerate(ZAs):
        print(f"\tRunning: Z={Z}, A={A} ({i_ZA+1} out of {len(ZAs)})")

        radionuclides["ZA"] = [Z, A]
        run_geometry(geometry, run_dict=radionuclides, run_type="radionuclides")




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

