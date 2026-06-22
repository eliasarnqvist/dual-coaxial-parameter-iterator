# Dual-coaxial-parameter-iterator

This repository contains code for simulating and evaluating dual-detector coincidence gamma-ray spectrometers. Contrary to the repository name, both generic coaxial and planar detectors are modeled. The purpose is to determine the detector minimum detectable activity (MDA) of detector setups. 

The MDA of a a pair of coincident gamma-rays is 

$$
MDA = \dfrac{L_D}{\varepsilon_{\gamma\gamma} I_{\gamma\gamma} t_{M}} ,
$$

where $L_D$ is the detection limit and is a function of the amount of background radiation, $\varepsilon_{\gamma\gamma}$ is the absolute full-energy peak gamma-ray coincidence detection efficiency, $I_{\gamma\gamma}$ is the coincidence gamma-ray emission intensity, and $t_{M}$ is the measurement time. Therefore, to caluclate the MDA, both background and efficiency need to be considered. In Geant4, radionuclide simulations are used to determine the product $\varepsilon_{\gamma\gamma} I_{\gamma\gamma}$ and background simulations are used to determine $L_D$. Python code is included for data analysis, which yields MDA as a function of a detector design parameter, facilitating parameteric design optimization of a detector concept. 

## Instructions

Clone the repository and (add micromamba instructions here)

### Geant4 simulations



### Data analysis


