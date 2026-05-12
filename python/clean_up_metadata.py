import json

metadata_file = "../geant4/output_background/metadata.json"
with open(metadata_file, "r") as f:
    metadata = json.load(f)

for key, value in metadata.items():
    metadata[key]["properties"]["SURE_background_total_flux"] = value["properties"]["SURE_background_total_flux"] / 100
    metadata[key]["properties"]["SURE_pseudo_time"] = value["properties"]["SURE_pseudo_time"] * 100

print(metadata)

# with open(metadata_file, "w") as f:
#     json.dump(metadata, f, indent=4)