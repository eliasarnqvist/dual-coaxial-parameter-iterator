def create_data_list(plotcard, metadata):
    # Make a "list mode" data storage scheme
    # A list with dictionaries
    # Every list entry is a specific data point with a specific detector design, radionuclide, and gamma ray(s)
    # The dictionary specifies which properties this data point has
    # Calculated quantities (like the MDAs) are added to this list later on

    data_list = []
    
    # First use the metadata to supply the geometry and radionuclide information
    for key, value in metadata.items():
        if value["type"] == "radionuclides":
            data_point = {
                "detector_type": value["properties"]["detector_type"],
                "detector_diameter": value["properties"]["detector_diameter"],
                "detector_length": value["properties"]["detector_length"],
                "detector_source_distance": value["properties"]["detector_source_distance"],
                "source_type": value["properties"]["source_type"],
                "ZA": [value["properties"]["Z"], value["properties"]["A"]],
            }
            if data_point not in data_list:
                data_list.append(data_point)
            else:
                # Somehow there is a duplicate!
                raise ValueError
        else:
            # Not metadata for a radionuclide
            pass
    
    # Next use the plotcard to supply the background file information
    new_data_list = []
    background_files = plotcard["background"]["background_file"]
    for data_point in data_list:
        for background_file in background_files:
            new_data_point = data_point.copy() # important with .copy() here!
            new_data_point["background_file"] = background_file
            new_data_list.append(new_data_point)
    data_list = new_data_list

    # Next use the plotcard to supply the gamma information
    new_data_list = []
    ZAs = plotcard["radionuclides"]["ZAs"]
    for data_point in data_list:
        for i, ZA in enumerate(ZAs):
            if data_point["ZA"] == ZA:
                # Check which energies are relevant for this radionuclide
                gammas_singles = plotcard["radionuclides"]["gammas_singles"][i]
                gammas_coincidences = plotcard["radionuclides"]["gammas_coincidences"][i]

                # Need a new data point for every single gamma ray energy
                for gamma_single in gammas_singles:
                    new_data_point = data_point.copy() # important with .copy() here!
                    new_data_point["analysis_type"] = "singles"
                    new_data_point["gamma"] = gamma_single
                    new_data_list.append(new_data_point)
                
                # Need a new data point for every pair of coincidence gamma rays
                for gamma_coincidence in gammas_coincidences:
                    new_data_point = data_point.copy() # important with .copy() here!
                    new_data_point["analysis_type"] = "coincidences"
                    new_data_point["gamma"] = gamma_coincidence
                    new_data_list.append(new_data_point)
            else:
                # This data point has a different radionuclide
                pass
    data_list = new_data_list

    # Next use the plotcard to supply the filter background information
    new_data_list = []
    filter_actives = plotcard["filter"]["filter_active"]
    for data_point in data_list:
        # Need a new data point for the filter activity options
        for filter_active in filter_actives:
            new_data_point = data_point.copy() # important with .copy() here!
            new_data_point["filter_active"] = filter_active
            new_data_list.append(new_data_point)
    data_list = new_data_list

    return data_list