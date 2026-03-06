import subprocess
import numpy as np
import time
import uuid
import json
import os

# %% Settings

# Choose what to do
use_filter_source = "true" # "true" or "false"

# The background files to iterate over
resource_folder = "resources/"
background_filenames = ["source_flux_avg_unshielded.csv", "source_flux_avg_shielded.csv"]
def calc_source_radius(diameter, length):
    R = np.sqrt((length+1.5+3+5+5+10+5+13.8/2)**2 + (diameter/2+2+3+1.5)**2)
    return R
# Detector parameters to iterate over
# detector_diameters = [50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0] # mm
detector_diameters = [85.0, 90.0] # mm
detector_lengths = "same as diameter"

# Number of runs per setting
runs_per_background = 5e8
# 1e7 = 1.1 min laptop 18t
# 5e7 = 5.4 min laptop 18t
# 1e8 = 10 min pelle 20t

build_folder = "build/"
macro_name = "autorun.mac"
output_folder = "output/"
number_of_threads = 32

# %% Run
print("Starting Geant4 autorun")
print("Number of threads: " + str(number_of_threads))
print("Settings: ")
print("\tdetector diameters: " + str(detector_diameters))
print("\tdetector lengths: " + str(detector_lengths))
print("\tbackground files: " + str(background_filenames))
print("\truns per background: " + str(int(runs_per_background)))

start_time = time.time()
time_interval = start_time

settings = [(bac, dia) 
            for bac in background_filenames
            for dia in detector_diameters]
            # for len in detector_lengths]

# Iterate over background
for i_b, (background_filename, diameter) in enumerate(settings):
    length = diameter

    print(f"Running setting {i_b + 1} of {len(settings)}: ")
    print(f"\tfile={background_filename}, diameter={diameter}, length={length}")

    flux_data = np.genfromtxt(resource_folder + background_filename, delimiter=",")
    E = flux_data[:, 0]
    flux = flux_data[:, 1:]
    background_total_flux = flux.sum()
    # Make cumulative distribution function
    flux_cdf = np.cumsum(flux) / background_total_flux
    # Save with all 9 decimals
    np.savetxt("resources/flux_cdf.dat", np.column_stack((E, flux_cdf)), fmt="%.9f\t%.9f")

    source_radius = calc_source_radius(diameter, length)
    
    pseudo_time = runs_per_background / (source_radius**2 * np.pi * background_total_flux)

    settings_name = "_dia_" + str(diameter) + "_len_" + str(length)
    settings_name += "_.root"

    # Make macro file first
    macro_content = ""
    macro_content += "/run/numberOfThreads " + str(number_of_threads) + "\n"
    file_name = output_folder + "threadoutput_" + str(i_b) + settings_name
    macro_content += "/E_file_settings/fileName " + file_name + "\n"
    macro_content += "/E_detector/detectorDiameter " + str(diameter) + "\n"
    macro_content += "/E_detector/detectorLength " + str(length) + "\n"
    macro_content += "/E_source/selectFilterSource " + use_filter_source + "\n"
    macro_content += "/run/reinitializeGeometry" + "\n"
    macro_content += "/run/initialize" + "\n"
    macro_content += "/E_source/sourceRadius " + str(source_radius) + "\n"
    macro_content += "/E_source/selectBackground true" + "\n"
    macro_content += "/E_source/selectFilterSource " + use_filter_source + "\n"
    macro_content += "/run/printProgress " + str(int(runs_per_background/10)) + "\n"
    macro_content += "/run/beamOn " + str(int(runs_per_background))

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
    process_root = "hadd -f " + output_file + " " + output_folder + "threadoutput_" + str(i_b) + "*.root"
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
        "type":"external_background_with_filter",
        "properties":{
            "model":"coaxial n-type v1",
            "background_file":background_filename,
            "diameter":diameter,
            "length":length,
            "sourcetype":"FOI filter v1",
            "runs":runs_per_background,
            "source_radius":source_radius,
            "background_total_flux":background_total_flux,
            "pseudo_time":pseudo_time,
            "threads":number_of_threads,
            "time":simulated_minutes,
            },
        }
    with open(output_folder + "metadata.json", "w") as f:
         json.dump(metadata, f, indent=4)

    print("\tDeleting temporary ROOT files...")
    process_delete = "rm " + output_folder + "threadoutput_" + str(i_b) + "*.root"
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


