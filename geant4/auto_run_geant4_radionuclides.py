import subprocess
import numpy as np
import time
import uuid
import json
import os

# %% Settings

# Choose what to do
use_filter_source = "true" # "true" or "false"

# The radionuclides to simulate (Z, A)
# ZAs = [(57, 140)]
# La-140, Ba-140, Mo-99, Te-132, I-132, I-131, I-133, I-135, Np-239
ZAs = [(57, 140), (56, 140), (42, 99), (52, 132), (53, 132), (53, 131), (53, 133), (53, 135), (93, 239)]
# Detector parameters to iterate over
detector_diameters = [50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0] # mm
detector_lengths = "same as diameter"

# Number of runs per setting
runs_per_radionuclide = 1e7
# 1e7 = 8 min per ZA pelle 20t

build_folder = "build/"
macro_name = "autorun.mac"
output_folder = "output/"
number_of_threads = 16

# %% Run
print("Starting Geant4 autorun")
print("Number of threads: " + str(number_of_threads))
print("Settings: ")
print("\tradionuclides: " + str(ZAs))
print("\tdetector diameters: " + str(detector_diameters))
print("\tdetector lengths: " + str(detector_lengths))
print("\truns per radionuclide: " + str(int(runs_per_radionuclide)))

start_time = time.time()
time_interval = start_time

settings = [(ZA, dia) 
            for ZA in ZAs
            for dia in detector_diameters]
            # for len in detector_lengths]

# Iterate over radionuclides
for i_s, (ZA, diameter) in enumerate(settings):
    # if i_s <= 2:
    #     continue
    
    length = diameter

    print(f"Running setting {i_s + 1} of {len(settings)}: ")
    print(f"\tZA={ZA}, diameter={diameter}, length={length}")

    Z = int(ZA[0])
    A = int(ZA[1])

    settings_name = "_Z_" + str(Z) + "_A_" + str(A)
    settings_name += "_dia_" + str(diameter) + "_len_" + str(length)
    settings_name += "_.root"

    # Make macro file first
    macro_content = ""
    macro_content += "/run/numberOfThreads " + str(number_of_threads) + "\n"
    file_name = output_folder + "threadoutput_" + str(i_s) + settings_name
    macro_content += "/E_file_settings/fileName " + file_name + "\n"
    macro_content += "/E_detector/detectorDiameter " + str(diameter) + "\n"
    macro_content += "/E_detector/detectorLength " + str(length) + "\n"
    macro_content += "/E_source/selectFilterSource " + use_filter_source + "\n"
    macro_content += "/run/reinitializeGeometry" + "\n"
    macro_content += "/run/initialize" + "\n"
    macro_content += "/process/had/rdm/thresholdForVeryLongDecayTime 1.0e+60 year" + "\n"
    macro_content += "/gun/particle ion" + "\n"
    macro_content += "/gun/ion " + str(Z) + " " + str(A) + " 0 0" + "\n"
    macro_content += "/process/had/rdm/nucleusLimits "+str(A)+" "+str(A)+" "+str(Z)+" "+str(Z)+"\n"
    macro_content += "/E_source/selectBackground false" + "\n"
    macro_content += "/E_source/selectFilterSource " + use_filter_source + "\n"
    macro_content += "/run/printProgress " + str(int(runs_per_radionuclide/10)) + "\n"
    macro_content += "/run/beamOn " + str(int(runs_per_radionuclide))

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
        "type":"radionuclides_in_filter",
        "properties":{
            "model":"coaxial n-type v1",
            "Z":Z,
            "A":A,
            "diameter":diameter,
            "length":length,
            "source_type":"FOI filter v1",
            "runs":runs_per_radionuclide,
            "threads":number_of_threads,
            "time":simulated_minutes,
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

print("All finished!")

end_time = time.time() # End timing
elapsed_minutes = (end_time - start_time) / 60
print(f"Total time spent: {elapsed_minutes:.2f} minutes")


