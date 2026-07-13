import uproot
import numpy as np
from analysis.analysis_helpers import ROI_analysis_1D
from analysis.analysis_helpers import ROI_analysis_2D
from analysis.analysis_helpers import calculate_effint
from analysis.analysis_helpers import calculate_B
from analysis.analysis_helpers import calculate_LD
from analysis.analysis_helpers import calculate_mda
from analysis.analysis_helpers import calculate_combined_mda
from analysis.analysis_helpers import calculate_B_filter


# TODO determine optimal measurement time for lowest mda
def analyze_files(data_list, plotcard, metadata, data_path):
    print("Starting analysis of data files!")

    # Quantities not needed inside a loop
    # TODO Optional determination of the measurement time for optimal mda
    measurement_time_hours = plotcard["analysis"]["measurement_time_hours"]

    # Need to iterate through the relevant files
    for i, (key, value) in enumerate(metadata.items()):

        # #s just for testing
        # if i > 4:
        #     continue

        run_type = value["type"]
        data_filename = value["filename"]
        data_filesize = value["file_size"] / (1024**2)
        print("\tOpening file " + str(i+1) + " out of " + str(len(metadata)) + ": " + str(data_filename) + " (" + run_type + f", {data_filesize:.2f} mb)...")

        data_filepath = data_path + "/" + data_filename
        data_file = uproot.open(data_filepath)
        # Single gamma rays
        tree_detector_a = data_file["E_a_list"]
        tree_detector_b = data_file["E_b_list"]
        events_single_a = tree_detector_a["energy_a"].array(library="np") * 1e3
        events_single_b = tree_detector_b["energy_b"].array(library="np") * 1e3
        # Coincidence gamma rays
        tree_detector_ab = data_file["E_ab_list"]
        events_coincidence_a = tree_detector_ab["energy_a"].array(library="np") * 1e3
        events_coincidence_b = tree_detector_ab["energy_b"].array(library="np") * 1e3

        # At the moment the data loading above is the bottleneck...
        # Cannot think of any real alternative to speed it up :/

        # Based on the file metadata, this should be the data point
        data_point_criteria = {
            "detector_type": value["properties"]["detector_type"],
            "detector_diameter": value["properties"]["detector_diameter"],
            "detector_length": value["properties"]["detector_length"],
            "detector_source_distance": value["properties"]["detector_source_distance"],
            "source_type": value["properties"]["source_type"],
        }

        if run_type == "radionuclides":
            # Based on the file, this should be the ZA of the data_point
            ZA_criteria = [value["properties"]["Z"], value["properties"]["A"]]
            data_point_criteria["ZA"] = ZA_criteria

            for j, ZA in enumerate(plotcard["radionuclides"]["ZAs"]):
                if ZA == ZA_criteria:
                    # Check which energies are relevant for this radionuclide
                    gammas_singles = plotcard["radionuclides"]["gammas_singles"][j]
                    gammas_coincidences = plotcard["radionuclides"]["gammas_coincidences"][j]

                    # Analyze single gamma rays
                    for E_gamma in gammas_singles:
                        if E_gamma == "combined":
                            # Handle this later
                            continue
                        else:
                            ROI_polynial_coeffs = plotcard["analysis"]["ROI_width_polynomial"]
                            ROI_width = np.polyval(ROI_polynial_coeffs, E_gamma)

                            counts_a = ROI_analysis_1D(events_single_a, E_gamma, ROI_width)
                            counts_b = ROI_analysis_1D(events_single_b, E_gamma, ROI_width)
                            print("\t\tZA: " + str(ZA) + ", gamma ray: " + str(E_gamma) + " keV, counts a: " + str(counts_a) + ", counts b: " + str(counts_b))

                            counts = counts_a + counts_b
                            events = value["properties"]["events"]

                            # Caluclate the effint
                            effint, effint_unc = calculate_effint(counts, events)

                        # Need to make a copy below so the keys for singles do not interfere with the keys for coincidences
                        this_data_point_criteria = data_point_criteria.copy()
                        this_data_point_criteria["analysis_type"] = "singles"
                        this_data_point_criteria["gamma"] = E_gamma

                        # Store the data in the data_list
                        for data_point in data_list:
                            if all(data_point[k] == v for k, v in this_data_point_criteria.items()):
                                # data_point["counts_radionuclide"] = counts
                                data_point["effint"] = [effint, effint_unc]

                    # Analyze coincidence gamma rays also
                    for Egamma1_Egamma2 in gammas_coincidences:
                        if Egamma1_Egamma2 == "combined":
                            # Handle this later
                            continue
                        else:
                            E_gamma1, E_gamma2 = Egamma1_Egamma2

                            ROI_polynial_coeffs = plotcard["analysis"]["ROI_width_polynomial"]
                            ROI_width1 = np.polyval(ROI_polynial_coeffs, E_gamma1)
                            ROI_width2 = np.polyval(ROI_polynial_coeffs, E_gamma2)

                            counts_a1b2 = ROI_analysis_2D(events_coincidence_a, events_coincidence_b, E_gamma1, E_gamma2, ROI_width1, ROI_width2)
                            counts_b1a2 = ROI_analysis_2D(events_coincidence_b, events_coincidence_a, E_gamma1, E_gamma2, ROI_width1, ROI_width2)
                            print("\t\tZA: " + str(ZA) + ", gamma ray 1: " + str(E_gamma1) + " keV, gamma ray 2: " + str(E_gamma2) + " keV , counts a1b2: " + str(counts_a1b2) + ", counts b2a1: " + str(counts_b1a2))

                            counts = counts_a1b2 + counts_b1a2
                            events = value["properties"]["events"]

                            # Caluclate the effint
                            effint, effint_unc = calculate_effint(counts, events)

                        # Need to make a copy below so the keys for singles do not interfere with the keys for coincidences
                        this_data_point_criteria = data_point_criteria.copy()
                        this_data_point_criteria["analysis_type"] = "coincidences"
                        this_data_point_criteria["gamma"] = Egamma1_Egamma2

                        # Store the data in the data_list
                        for data_point in data_list:
                            if all(data_point[k] == v for k, v in this_data_point_criteria.items()):
                                # data_point["counts_radionuclide"] = counts
                                data_point["effint"] = [effint, effint_unc]

        elif run_type == "background":
            # Based on the file, this should be the background file of the data point
            background_file_criteria = value["properties"]["background_file"]
            data_point_criteria["background_file"] = background_file_criteria

            # Iterate through all radionuclides and all gamma rays
            for j, ZA in enumerate(plotcard["radionuclides"]["ZAs"]):
                print("\t\tRadionuclide " + str(j+1) + " out of " + str(len(plotcard["radionuclides"]["ZAs"])) + f": Z={ZA[0]}, A={ZA[1]}")

                data_point_criteria["ZA"] = ZA
                
                # Check which energies are relevant for this radionuclide
                gammas_singles = plotcard["radionuclides"]["gammas_singles"][j]
                gammas_coincidences = plotcard["radionuclides"]["gammas_coincidences"][j]

                # Analyze single gamma rays
                for E_gamma in gammas_singles:
                    if E_gamma == "combined":
                        # Handle this later
                        continue
                    else:
                        ROI_polynial_coeffs = plotcard["analysis"]["ROI_width_polynomial"]
                        ROI_width = np.polyval(ROI_polynial_coeffs, E_gamma)

                        counts_a = ROI_analysis_1D(events_single_a, E_gamma, ROI_width)
                        counts_b = ROI_analysis_1D(events_single_b, E_gamma, ROI_width)
                        print("\t\t\tGamma ray: " + str(E_gamma) + " keV, counts a: " + str(counts_a) + ", counts b: " + str(counts_b))

                        counts = counts_a + counts_b
                        events = value["properties"]["events"]
                        pseudo_time = value["properties"]["SURE_pseudo_time"]

                        # Calculate the number of background counts
                        B, B_unc = calculate_B(counts, events, pseudo_time, measurement_time_hours)
                        
                    # Need to make a copy below so the keys for singles do not interfere with the keys for coincidences
                    this_data_point_criteria = data_point_criteria.copy()
                    this_data_point_criteria["analysis_type"] = "singles"
                    this_data_point_criteria["gamma"] = E_gamma

                    # Store the data in the data_list
                    for data_point in data_list:
                        if all(data_point[k] == v for k, v in this_data_point_criteria.items()):
                            # data_point["counts_radionuclide"] = counts
                            data_point["B"] = [B, B_unc]
                    
                # Analyze coincidence gamma rays also
                for Egamma1_Egamma2 in gammas_coincidences:
                    if Egamma1_Egamma2 == "combined":
                        # Handle this later
                        continue
                    else:
                        E_gamma1, E_gamma2 = Egamma1_Egamma2

                        ROI_polynial_coeffs = plotcard["analysis"]["ROI_width_polynomial"]
                        ROI_width1 = np.polyval(ROI_polynial_coeffs, E_gamma1)
                        ROI_width2 = np.polyval(ROI_polynial_coeffs, E_gamma2)

                        counts_a1b2 = ROI_analysis_2D(events_coincidence_a, events_coincidence_b, E_gamma1, E_gamma2, ROI_width1, ROI_width2)
                        counts_b1a2 = ROI_analysis_2D(events_coincidence_b, events_coincidence_a, E_gamma1, E_gamma2, ROI_width1, ROI_width2)
                        print("\t\t\tGamma ray 1: " + str(E_gamma1) + " keV, gamma ray 2: " + str(E_gamma2) + " keV , counts a1b2: " + str(counts_a1b2) + ", counts b2a1: " + str(counts_b1a2))

                        counts = counts_a1b2 + counts_b1a2
                        events = value["properties"]["events"]
                        pseudo_time = value["properties"]["SURE_pseudo_time"]

                        # Calculate the number of background counts
                        B, B_unc = calculate_B(counts, events, pseudo_time, measurement_time_hours)

                    # Need to make a copy below so the keys for singles do not interfere with the keys for coincidences
                    this_data_point_criteria = data_point_criteria.copy()
                    this_data_point_criteria["analysis_type"] = "coincidences"
                    this_data_point_criteria["gamma"] = Egamma1_Egamma2

                    # Store the data in the data_list
                    for data_point in data_list:
                        if all(data_point[k] == v for k, v in this_data_point_criteria.items()):
                            # data_point["counts_radionuclide"] = counts
                            data_point["B"] = [B, B_unc]

        elif run_type == "filter":
            # A criteria of applicability is as follows
            data_point_criteria["filter_active"] = True
            # Identify the ZA of the simulation
            filter_ZA = [value["properties"]["Z"], value["properties"]["A"]]

            # Go through the radionuclides and analyze the filter contribution
            for j, ZA in enumerate(plotcard["radionuclides"]["ZAs"]):
                print("\t\tRadionuclide " + str(j+1) + " out of " + str(len(plotcard["radionuclides"]["ZAs"])) + f": Z={ZA[0]}, A={ZA[1]}")

                data_point_criteria["ZA"] = ZA

                # Check which energies are relevant for this radionuclide
                gammas_singles = plotcard["radionuclides"]["gammas_singles"][j]
                gammas_coincidences = plotcard["radionuclides"]["gammas_coincidences"][j]

                # Analyze single gamma rays
                for E_gamma in gammas_singles:
                    if E_gamma == "combined":
                        # Handle this later
                        continue
                    else:
                        ROI_polynial_coeffs = plotcard["analysis"]["ROI_width_polynomial"]
                        ROI_width = np.polyval(ROI_polynial_coeffs, E_gamma)

                        counts_a = ROI_analysis_1D(events_single_a, E_gamma, ROI_width)
                        counts_b = ROI_analysis_1D(events_single_b, E_gamma, ROI_width)
                        print("\t\t\tGamma ray: " + str(E_gamma) + " keV, counts a: " + str(counts_a) + ", counts b: " + str(counts_b))

                        counts = counts_a + counts_b
                        events = value["properties"]["events"]
                        ZA_Bq_filter = next((x for x in plotcard["filter"]["ZAs_Bq_filter"] if x[:2] == filter_ZA), None)
                        Bq = ZA_Bq_filter[2]
                        Bq_unc = ZA_Bq_filter[3]

                        # Calculate the number of background counts
                        B_filter, B_filter_unc = calculate_B_filter(counts, events, [Bq, Bq_unc], measurement_time_hours)
                    
                    # Need to make a copy below so the keys for singles do not interfere with the keys for coincidences
                    this_data_point_criteria = data_point_criteria.copy()
                    this_data_point_criteria["analysis_type"] = "singles"
                    this_data_point_criteria["gamma"] = E_gamma
                    
                    # Store the data in the data_list
                    for data_point in data_list:
                        if all(data_point[k] == v for k, v in this_data_point_criteria.items()):
                            if "B_filter" in data_point:
                                new_B_filter = data_point["B_filter"][0] + B_filter
                                new_B_filter_unc = np.sqrt(np.power(data_point["B_filter"][1], 2) + B_filter_unc)
                                data_point["B_filter"] = [new_B_filter, new_B_filter_unc]
                            else:
                                data_point["B_filter"] = [B_filter, B_filter_unc]
                
                # Analyze single gamma rays
                # Analyze coincidence gamma rays also
                for Egamma1_Egamma2 in gammas_coincidences:
                    if Egamma1_Egamma2 == "combined":
                        # Handle this later
                        continue
                    else:
                        E_gamma1, E_gamma2 = Egamma1_Egamma2

                        ROI_polynial_coeffs = plotcard["analysis"]["ROI_width_polynomial"]
                        ROI_width1 = np.polyval(ROI_polynial_coeffs, E_gamma1)
                        ROI_width2 = np.polyval(ROI_polynial_coeffs, E_gamma2)

                        counts_a1b2 = ROI_analysis_2D(events_coincidence_a, events_coincidence_b, E_gamma1, E_gamma2, ROI_width1, ROI_width2)
                        counts_b1a2 = ROI_analysis_2D(events_coincidence_b, events_coincidence_a, E_gamma1, E_gamma2, ROI_width1, ROI_width2)
                        print("\t\t\tGamma ray 1: " + str(E_gamma1) + " keV, gamma ray 2: " + str(E_gamma2) + " keV , counts a1b2: " + str(counts_a1b2) + ", counts b2a1: " + str(counts_b1a2))

                        counts = counts_a1b2 + counts_b1a2
                        events = value["properties"]["events"]
                        ZA_Bq_filter = next((x for x in plotcard["filter"]["ZAs_Bq_filter"] if x[:2] == filter_ZA), None)
                        Bq = ZA_Bq_filter[2]
                        Bq_unc = ZA_Bq_filter[3]

                        # Calculate the number of background counts
                        B_filter, B_filter_unc = calculate_B_filter(counts, events, [Bq, Bq_unc], measurement_time_hours)
                    
                    # Need to make a copy below so the keys for singles do not interfere with the keys for coincidences
                    this_data_point_criteria = data_point_criteria.copy()
                    this_data_point_criteria["analysis_type"] = "coincidences"
                    this_data_point_criteria["gamma"] = Egamma1_Egamma2
                    
                    # Store the data in the data_list
                    for data_point in data_list:
                        if all(data_point[k] == v for k, v in this_data_point_criteria.items()):
                            if "B_filter" in data_point:
                                new_B_filter = data_point["B_filter"][0] + B_filter
                                new_B_filter_unc = np.sqrt(np.power(data_point["B_filter"][1], 2) + B_filter_unc)
                                data_point["B_filter"] = [new_B_filter, new_B_filter_unc]
                            else:
                                data_point["B_filter"] = [B_filter, B_filter_unc]
        else:
            pass

    # Now we can calculate the detection limit for the data points
    for data_point in data_list:
        if data_point["gamma"] == "combined":
                # Handle this later
                continue

        # Check if background and/or filter background is present for this data point
        B_tot = 0
        B_tot_unc = 0
        if ("B" in data_point):
            B_tot += data_point["B"][0]
            B_tot_unc = np.sqrt(np.power(B_tot_unc, 2) + np.power(data_point["B"][1], 2))
        if ("B_filter" in data_point):
            B_tot += data_point["B_filter"][0]
            B_tot_unc = np.sqrt(np.power(B_tot_unc, 2) + np.power(data_point["B_filter"][1], 2))

        LD, LD_unc = calculate_LD([B_tot, B_tot_unc])

        data_point["LD"] = [LD, LD_unc]

    # Now it is possible to caluclate the MDA for each data point
    for data_point in data_list:
        # Check that both efficiency and intensity information is present for this data point
        if ("effint" in data_point) and ("LD" in data_point):
            # Get the effint and LD information for this data point
            effint, effint_unc = data_point["effint"]
            LD, LD_unc = data_point["LD"]

            if data_point["gamma"] == "combined":
                # Handle this later
                continue
            
            mda, mda_unc = calculate_mda([LD, LD_unc], [effint, effint_unc], measurement_time_hours, t12=0)

            data_point["mda"] = [mda, mda_unc]
        else:
            # This is the case when either of effint or LD are not calculated
            data_point["mda"] = [1000.0, 1000.0]
    
    # Now calculate the combined MDA entries where applicable
    for data_point in data_list:
        # Check if this data point is for combined signatures
        if data_point["gamma"] != "combined":
            continue

        mda_list = []

        # Check which entries have the same properties as this data point
        data_point_criteria = {
            "detector_type": data_point["detector_type"],
            "detector_diameter": data_point["detector_diameter"],
            "detector_length": data_point["detector_length"],
            "detector_source_distance": data_point["detector_source_distance"],
            "source_type": data_point["source_type"],
            "ZA": data_point["ZA"],
            "background_file": data_point["background_file"],
            "analysis_type": data_point["analysis_type"],
            "filter_active": data_point["filter_active"],
        }
        for data_point_candidate in data_list:
            if all(data_point_candidate[k] == v for k, v in data_point_criteria.items()) and data_point_candidate["gamma"] != "combined":
                mda, mda_unc = data_point_candidate["mda"]
                mda_list.append([mda, mda_unc])

        combined_mda, combined_mda_unc = calculate_combined_mda(mda_list)
        data_point["mda"] = [combined_mda, combined_mda_unc]
    
    return data_list

