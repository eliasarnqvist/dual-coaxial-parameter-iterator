import json

# metadata_file = "output_pelle/metadata.json"
metadata_file = "pelle2/output/metadata.json"

with open(metadata_file, "r") as f:
    metadata = json.load(f)

for key, value in metadata.items():

    # metadata[key]["properties"]["SURE_background_total_flux"] = value["properties"]["SURE_background_total_flux"] / 100
    # metadata[key]["properties"]["SURE_pseudo_time"] = value["properties"]["SURE_pseudo_time"] * 100

    # Fix source type
    try:
        select_filter_source = metadata[key]["properties"]["select_filter_source"]
        if select_filter_source == True:
            metadata[key]["properties"]["source_type"] = 1
        else:
            metadata[key]["properties"]["source_type"] = 0
        metadata[key]["properties"].pop("select_filter_source")
    except:
        pass

    # Fix SURE radius
    try:
        if metadata[key]["type"] == "radionuclides":
            metadata[key]["properties"].pop("source_SURE_radius")
    except:
        pass

    # Fix model or detector type
    try:
        if metadata[key]["properties"]["select_ntype_instead_of_ptype"] == True:
            metadata[key]["properties"]["detector_type"] = 1
        else:
            metadata[key]["properties"]["detector_type"] = 0
        metadata[key]["properties"].pop("model")
        metadata[key]["properties"].pop("select_ntype_instead_of_ptype")
    except:
        pass

    # Background
    try:
        if metadata[key]["type"] == "external_background":
            metadata[key]["type"] = "background"
        metadata[key]["properties"].pop("SURE_background_total_flux")
    except:
        pass

    try:
        if metadata[key]["type"] == "background":
            metadata[key]["properties"]["SURE_radius"] = metadata[key]["properties"]["source_SURE_radius"]
            metadata[key]["properties"].pop("source_SURE_radius")
    except:
        pass

    try:
        metadata[key]["properties"]["background_file"] = "resources/" + metadata[key]["properties"]["background_file"]
    except:
        pass

    try:
        metadata[key]["properties"]["time_minutes"] = metadata[key]["properties"]["time"]
        metadata[key]["properties"].pop("time")
    except:
        pass

    # Filter
    try:
        if metadata[key]["type"] == "radionuclides":
            metadata[key]["type"] = "filter"
    except:
        pass


print(metadata)

# with open(metadata_file, "w") as f:
#     json.dump(metadata, f, indent=4)