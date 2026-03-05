# Dual-coaxial-parameter-iterator

Coincidence gamma-ray spectrometry can achieve greatly improved minimum detectable activities (MDA), compared to many conventional detector systems. However, this relies on that the coincidence detector system has been carefully designed. This repository contains a Geant4 application that models a coincidence gamma-ray spectrometer based on two generic coaxial HPGe detectors. The design parameters of the detectors can then be varied, as specified by the user. The goal is to determine which design parameters result in the best performance. 

The MDA of a a pair of coincident gamma-rays is 

$$
MDA = \dfrac{L_D}{\varepsilon_{\gamma\gamma} \, I_{\gamma\gamma} \, t} ,
$$

where $L_D$ is the detection limit and is a function of the amount of background radiation, $\varepsilon_{\gamma\gamma}$ is the absolute full-energy peak gamma-ray coincidence detection efficiency, $I_{\gamma\gamma}$ is the coincidence gamma-ray emission intensity, and $t$ is the measurement time. To caluclate the MDA, two simulations are performed. A radionuclide simulation is used to determine $\varepsilon_{\gamma\gamma} \, I_{\gamma\gamma}$ and a background simulation is used to determine $L_D$. 

Python code is also included for data analysis, which yields MDA as a function of a detector design parameter. In other words, facilitating parameteric design optimization of a detector concept. 

## Instructions

Will come soon. 

## Results

Will come soon. 

## Limitations

Will come soon. 
