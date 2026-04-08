import subprocess
import numpy as np
import time
import uuid
import json
import os

# %% Settings

# The background files to iterate over
resource_folder = "resources/"
background_filenames = ["source_flux_avg_unshielded.csv"]
# For calculating the SURE model radius
def calc_source_radius(diameter, length, select_filter_source):
    if select_filter_source:
        source_thickness = 13.8/2
    else:
        source_thickness = 0
    R = np.sqrt((length+1.5+4+5+5+5+5+source_thickness)**2 + (diameter/2+2+4+1.5)**2)
    return R

# Detector parameters to iterate over (in mm)
detector_diameters = [40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0]
detector_lengths = [-1] # if -1 then length=diameter
source_distances = [0]
select_n_type_instead_of_p_type = False
select_filter_source = True
source_type = 2 # 0=point, 1=filter, 2=SURE
source_SURE_radius = None

# Number of events per run
events_per_background = 1e9
# 1e8 = 10 min pelle 20t

build_folder = "build/"
macro_name = "autorun.mac"
output_folder = "output/"
number_of_threads = 32

# %% Run
print("Starting Geant4 autorun")
print("Number of threads: " + str(number_of_threads))
print("Settings: ")
print("\tbackground files: " + str(background_filenames))
print("\tdetector diameters: " + str(detector_diameters))
print("\tdetector lengths: " + str(detector_lengths))
print("\tsource distances: " + str(source_distances))
print("\tn-type (True) or p-type (False): " + str(select_n_type_instead_of_p_type))
print("\tfilter source (True) or point source (False): " + str(select_filter_source))
print("\truns per background: " + str(int(events_per_background)))

start_time = time.time()
time_interval = start_time

settings = [(bac, dia, len, dis) 
            for bac in background_filenames
            for dia in detector_diameters
            for len in detector_lengths
            for dis in source_distances]

# Iterate over background
for i_s, (background_filename, diameter, length, distance) in enumerate(settings):
    # Optional: skip some iterations
    # if i_s <= 2:
    #     continue
    
    if length == -1:
        length = diameter
    
    print(f"Running setting {i_s + 1} of {len(settings)}: ")
    print(f"\tfile={background_filename}, diameter={diameter}, length={length}")

    # Background source spectrum for SURE model
    flux_data = np.genfromtxt(resource_folder + background_filename, delimiter=",")
    E = flux_data[:, 0]
    flux = flux_data[:, 1:]
    background_total_flux = flux.sum()
    # Make cumulative distribution function
    flux_cdf = np.cumsum(flux) / background_total_flux
    # Save with all 9 decimals
    np.savetxt("resources/flux_cdf.dat", np.column_stack((E, flux_cdf)), fmt="%.9f\t%.9f")

    if source_type == 2:
        source_SURE_radius = calc_source_radius(diameter, length, select_filter_source)
    
    pseudo_time = events_per_background / (source_SURE_radius**2 * np.pi * background_total_flux)

    settings_name = "_background_"
    settings_name += "_dia_" + str(diameter) + "_len_" + str(length) + "_dis_" + str(distance)
    settings_name += "_.root"

    # Make macro file first
    macro_content = ""
    macro_content += "/run/numberOfThreads " + str(number_of_threads) + "\n"

    file_name = output_folder + "threadoutput_" + str(i_s) + settings_name
    macro_content += "/E_file_settings/fileName " + file_name + "\n"

    macro_content += "/E_detector/detectorDiameter " + str(diameter) + "\n"
    macro_content += "/E_detector/detectorLength " + str(length) + "\n"
    macro_content += "/E_detector/sourceDistance " + str(distance) + "\n"
    macro_content += "/E_detector/selectNTypeInsteadOfPType " + str(select_n_type_instead_of_p_type) + "\n"
    macro_content += "/E_detector/selectFilterSource " + str(select_filter_source) + "\n"

    macro_content += "/run/reinitializeGeometry" + "\n"
    macro_content += "/run/initialize" + "\n"
    macro_content += "/process/had/rdm/thresholdForVeryLongDecayTime 1.0e+60 year" + "\n"

    macro_content += "/E_source/sourceType " + str(source_type) + "\n"
    macro_content += "/E_source/sourceRadiusSURE " + str(source_SURE_radius) + "\n"

    macro_content += "/run/printProgress " + str(int(events_per_background/10)) + "\n"
    macro_content += "/run/beamOn " + str(int(events_per_background))

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
        "type":"external_background",
        "properties":{
            "model":"coaxial_v2",
            "background_file":background_filename,
            "detector_diameter":diameter,
            "detector_length":length,
            "source_distance":distance,
            "select_n_type_instead_of_ptype":select_n_type_instead_of_p_type,
            "select_filter_source":select_filter_source,
            "source_type":source_type,
            "source_SURE_radius":source_SURE_radius,
            "SURE_background_total_flux":background_total_flux,
            "SURE_pseudo_time":pseudo_time,
            "source_type":"FOI_filter_v1",
            "events":events_per_background,
            "threads":number_of_threads,
            "time":simulated_minutes,
            "throughput":events_per_background/(simulated_minutes*60*number_of_threads),
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
