# Cut the metadata according to the cuts in the plotcard
def cut_metadata(plotcard, metadata):
    # Use the selection cuts in the plotcard to remove irrelevant metadata
    for key in list(metadata):
        value = metadata[key]
        
        keep_this_key = True

        # Check that the geometry matches (for all runs)
        for geometry_setting in ["detector_type", "detector_diameter", "detector_length", "detector_source_distance", "source_type"]:
            if len(plotcard["geometry"][geometry_setting]) == 1 and plotcard["geometry"][geometry_setting][0] == "any":
                pass
            elif value["properties"][geometry_setting] not in plotcard["geometry"][geometry_setting]:
                keep_this_key = False
            else:
                pass

        # Check that the radionuclide/background settings match (only for applicable runs)
        if value["type"] == "radionuclides":
            ZA = [value["properties"]["Z"], value["properties"]["A"]]
            if ZA not in plotcard["radionuclides"]["ZAs"]:
                keep_this_key = False
        elif value["type"] == "background":
            if value["properties"]["background_file"] not in plotcard["background"]["background_file"]:
                keep_this_key = False
        elif value["type"] == "filter":
            ZA = [value["properties"]["Z"], value["properties"]["A"]]
            ZAs_filter = [plotcard["filter"]["ZAs_Bq_filter"][i][0:2] for i in range(len(plotcard["filter"]["ZAs_Bq_filter"]))]
            filter_active = plotcard["filter"]["filter_active"]
            if (ZA not in ZAs_filter) or (True not in filter_active):
                keep_this_key = False
        else:
            raise ValueError
        
        # Remove the key from the metadata dictionary if it is not interesting
        if keep_this_key == False:
            metadata.pop(key)
    
    # Return the shortened metadata dictionary
    return metadata
