# import numpy as np
# from uncertainties import ufloat
# import uproot
import argparse
# from numba import jit
import json

# Parser for adding arguments
parser = argparse.ArgumentParser(prog="DetermineEffInt",
                                 description="Analyze .root files to determine the product of efficiency and intensity",
                                 epilog="Elias Arnqvist, 2026, Uppsala University",
                                 add_help=True)
parser.add_argument("-Z", type=int, required=True, help="Atomic number (Z) of radionuclide")
parser.add_argument("-A", type=int, required=True, help="Mass number (A) of radionuclide")
parser.add_argument("-d", "--data", type=str, required=False, default="../../geant4/output_pelle/", help="Path to directory with root files")
parser.add_argument("-m", "--metadata", type=str, required=False, default="../../geant4/output_pelle/metadata.json", help="Path to metadata file")
parser.add_argument("-g", "--gammas", type=str, required=False, default="gammas.json", help="Path to gamma ray information")
parser.add_argument("-s", "--save", type=str, required=False, default="/efficiency_intensity.json", help="Path to save output data")
args = parser.parse_args()

# Open gamma ray file
gammas_filepath = args.gammas
print("Opening gamma file: " + gammas_filepath)
with open(gammas_filepath, "r") as f:
    gammas = json.load(f)

# Open metadata filels
metadata_filepath = args.metadata
print("Opening metadata file: " + metadata_filepath)
with open(metadata_filepath, "r") as f:
    metadata = json.load(f)

# Look




# parser.add_argument("-s", "--save", type=str, default="figure", help="save path (default: 'figure')")
# parser.add_argument("-imin", "--imin", type=int, default=0.0, help="Minimum gamma intensity (default: 0.0)")
# parser.add_argument("-imax", "--imax", type=int, default=100.0, help="Maximum gamma intensity (default: 100.0)")
# args = parser.parse_args()


# @jit(nopython=True)
# def ROI_2D_square (events, Eg1, Eg2, ROI):

#     # ROI is the full width of the ROI
#     # Half of the ROI is the distance from the gamma peak to the ROI edge
#     dE = ROI / 2

#     branch_a = "energy_a"
#     branch_b = "energy_b"

#     Ea = events[branch_a]
#     Eb = events[branch_b]

#     cond_a1 = np.logical_and(Ea > Eg1-dE, Ea < Eg1+dE)
#     cond_b2 = np.logical_and(Eb > Eg2-dE, Eb < Eg2+dE)
#     cond_coincidence_a1b2 = np.logical_and(cond_a1, cond_b2)

#     cond_a2 = np.logical_and(Ea > Eg2-dE, Ea < Eg2+dE)
#     cond_b1 = np.logical_and(Eb > Eg1-dE, Eb < Eg1+dE)
#     cond_coincidence_a2b1 = np.logical_and(cond_a2, cond_b1)

#     return

# def ROI_1D (events, Eg, ROI):

#     dE = ROI / 2

