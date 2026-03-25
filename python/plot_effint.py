import matplotlib.pyplot as plt
import json


metadata_peakfinfo_filepath = "metadata_peakinfo.json"

with open(metadata_peakfinfo_filepath, "r") as f:
    metadata_peakfinfo = json.load(f)

ZAs_metadata = [(value["properties"]["Z"], value["properties"]["A"]) for value in metadata_peakfinfo.values() if value["type"] == "radionuclides"]
ZAs_metadata = list(set(ZAs_metadata))


plt.close('all')
inch_to_mm = 25.4

for i, (Z, A) in enumerate(ZAs_metadata):
    
    x = []
    y = []

    for i, (key, value) in enumerate(metadata_peakfinfo.items()):
            if value["type"] != "radionuclides":
                continue
            if value["properties"]["Z"] != Z or value["properties"]["A"] != A:
                continue
            if value["properties"]["select_n_type_instead_of_ptype"] != True:
                continue

            x.append(value["properties"]["detector_diameter"])


            y.append(value["peakinfo"]["detector_diameter"])


    
    fig, ax = plt.subplots(1, 1, figsize=(100/inch_to_mm,100/inch_to_mm))









