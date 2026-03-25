import argparse
import matplotlib.pyplot as plt
import json
import numpy as np


def plot_singles(ZAs):
    for i, (Z, A) in enumerate(ZAs):
        fig, ax = plt.subplots(1, 1, figsize=(100/inch_to_mm,100/inch_to_mm))

        x = []
        y = []
        y_unc = []

        for i, (key, value) in enumerate(metadata_peakfinfo.items()):
            if value["type"] != "radionuclides":
                continue
            if value["properties"]["Z"] != Z or value["properties"]["A"] != A:
                continue
            # if value["properties"]["select_n_type_instead_of_ptype"] != True:
            #     continue

            counts = value["peakdata"]["singles"]["counts_a"][1]

            events = value["properties"]["runs"]

            effint = counts / events
            effint_unc = np.sqrt(counts) / events

            x.append(value["properties"]["detector_diameter"])


            y.append(effint)
            y_unc.append(effint_unc)

        ax.errorbar(x, y, yerr=y_unc, fmt=".", ls="-", markersize=4, lw=1, capsize=2, capthick=1)
        # plt.show(block=False)

        plt.tight_layout(pad = 0.2)
        fig.subplots_adjust(hspace=0, wspace=0)
        plt.savefig("figures/test.jpg", dpi=300)


    return


def plot_coincidences():
    return


# Parser for adding arguments
parser = argparse.ArgumentParser(prog="plot_effint",
                                 description="Plot the product of efficiency and intensity",
                                 epilog="Elias Arnqvist, 2026, Uppsala University",
                                 add_help=True)
# parser.add_argument("-x", type=str, required=True, default="None", help="What to plot on the x-axis")
# parser.add_argument("-y", type=str, required=True, default="None", help="What to plot on the y-axis")
parser.add_argument("-Z", type=int, required=True, default=0, help="Atomic number (Z) of radionuclide")
parser.add_argument("-A", type=int, required=True, default=0, help="Mass number (A) of radionuclide")
parser.add_argument("-mp", "--metadata_peakinfo", type=str, required=False, default="metadata_peakinfo.json", help="Path to metadata and peakinfo file")
args = parser.parse_args()

metadata_peakfinfo_filepath = args.metadata_peakinfo
with open(metadata_peakfinfo_filepath, "r") as f:
    metadata_peakfinfo = json.load(f)

ZAs_metadata = [(value["properties"]["Z"], value["properties"]["A"]) for value in metadata_peakfinfo.values() if value["type"] == "radionuclides"]
ZAs_metadata = list(set(ZAs_metadata))

Z = args.Z
A = args.A
print("Using specified ZA: Z=" + str(Z) + ", A=" + str(A))
if (Z, A) not in ZAs_metadata:
    raise Exception("Used Z and A not found in metadata")

plt.close('all')
inch_to_mm = 25.4

plot_singles([(Z, A)])




