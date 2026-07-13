import numpy as np


def ROI_analysis_1D(event_data, E_gamma, ROI_width, events, ROI_background_extension_factor=1):
    # Energy distance is half of the ROI size
    dE = (ROI_width/2) * ROI_background_extension_factor

    # Count peak counts
    cond = np.logical_and(event_data > E_gamma-dE, event_data < E_gamma+dE)
    counts = cond.sum()
    counts_unc = np.sqrt(counts * (1 - counts/events))

    # Correct for ROI_background_extension_factor
    counts /= np.power(ROI_background_extension_factor, 2)
    counts_unc /= np.power(ROI_background_extension_factor, 2)

    return int(counts), float(counts_unc)


def ROI_analysis_2D(event_data_a, event_data_b, E_gamma1, E_gamma2, ROI_width1, ROI_width2, events, ROI_background_extension_factor=1):
    # Energy distance is half of the ROI size
    dE_1 = (ROI_width1/2) * ROI_background_extension_factor
    dE_2 = (ROI_width2/2) * ROI_background_extension_factor

    # Count peak counts
    cond_1 = np.logical_and(event_data_a > E_gamma1-dE_1, event_data_a < E_gamma1+dE_1)
    cond_2 = np.logical_and(event_data_b > E_gamma2-dE_2, event_data_b < E_gamma2+dE_2)
    cond = np.logical_and(cond_1, cond_2)
    counts = cond.sum()
    counts_unc = np.sqrt(counts * (1 - counts/events))

    # Correct for ROI_background_extension_factor
    counts /= np.power(ROI_background_extension_factor, 2)
    counts_unc /= np.power(ROI_background_extension_factor, 2)

    return int(counts), float(counts_unc)


def calculate_effint(counts_, events):
    # Counts and uncertainty of counts
    counts, counts_unc = counts_

    # Calculate efficiency and intensity
    effint = counts / events

    # Use binomial instead of poisson approximation
    effint_unc = counts_unc / events

    return float(effint), float(effint_unc)


def calculate_B(counts_, events, pseudo_time, measurement_time_hours):
    # Counts and uncertainty of counts
    counts, counts_unc = counts_

    # Convert to seconds
    measurement_time = measurement_time_hours * 3600
    # Background count rate
    b = counts / pseudo_time
    b_unc = counts_unc / pseudo_time
    # Background counts during measurement time
    B = b * measurement_time
    B_unc = b_unc * measurement_time
    return float(B), float(B_unc)


def calculate_B_filter(counts_, events, Bq_, measurement_time_hours):
    # Activity of radionuclide in filter
    Bq, Bq_unc = Bq_
    # Counts and uncertainty of counts
    counts, counts_unc = counts_

    # Convert to seconds
    measurement_time = measurement_time_hours * 3600
    # Calculate efficiency and intensity
    effint = counts / events
    # Use binomial instead of poisson approximation
    effint_unc = counts_unc / events
    # Background counts
    B_filter = effint * Bq * measurement_time
    B_filter_unc = np.sqrt(np.power(Bq*measurement_time*effint_unc, 2) + np.power(effint*measurement_time*Bq_unc, 2))

    return float(B_filter), float(B_filter_unc)


def calculate_LD(B_):
    # Background counts
    B, B_unc = B_

    if B == 0:
        print("Encountered 0 background!!!")
        raise ValueError
    # Detection limit according to Currie
    LD = 2.71 + 4.65*np.sqrt(B)
    # Calculate the uncertainty
    LD_unc = 4.65 * (0.5/np.sqrt(B)) * B_unc

    return float(LD), float(LD_unc)


def calculate_mda(LD_, effint_, measurement_time_hours, t12=0):
    # Convert to seconds
    tM = measurement_time_hours * 3600
    # Extract uncertainties
    LD, LD_unc = LD_
    effint, effint_unc = effint_
    
    # Calculate the minimum detectable activity
    mda = LD / (effint * tM)
    mda_unc = np.sqrt(np.power((LD_unc)/(effint*tM), 2) + np.power((LD*effint_unc)/(np.power(effint, 2)*tM), 2))

    # Optional correction for decay during measurement
    if t12 != 0:
        lambdaa = np.log(2)/t12 # extra a to respect Python reserved word
        decay_correction = (lambdaa * tM) / (1 - np.exp(-lambdaa * tM))
        mda *= decay_correction
        mda_unc *= decay_correction

    return float(mda), float(mda_unc)


def calculate_combined_mda(mda_list):

    # Combined MDA added in inverse quadrature
    sum_1 = 0
    sum_2 = 0
    for mda, mda_unc in mda_list:
        sum_1 += 1 / np.power(mda, 2)
        sum_2 += np.power(mda_unc, 2) / np.power(mda, 6)
    
    mda_combined = np.sqrt(1 / sum_1)

    mda_combined_unc = np.power(mda_combined, 3) * np.sqrt(sum_2)

    return float(mda_combined), float(mda_combined_unc)

