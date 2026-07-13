import numpy as np


# TODO ROI averaging: inflate by 4x and then divide by 16?
def ROI_analysis_1D(events, E_gamma, ROI_width):
    # Energy distance is half of the ROI size
    dE = ROI_width/2

    # Count peak counts
    cond = np.logical_and(events > E_gamma-dE, events < E_gamma+dE)
    counts = cond.sum()
    return int(counts)


# TODO ROI averaging: inflate by 4x and then divide by 16?
def ROI_analysis_2D(events_a, events_b, E_gamma1, E_gamma2, ROI_width1, ROI_width2):
    # Energy distance is half of the ROI size
    dE_1 = ROI_width1/2
    dE_2 = ROI_width2/2

    # Count peak counts
    cond_1 = np.logical_and(events_a > E_gamma1-dE_1, events_a < E_gamma1+dE_1)
    cond_2 = np.logical_and(events_b > E_gamma2-dE_2, events_b < E_gamma2+dE_2)
    cond = np.logical_and(cond_1, cond_2)
    counts = cond.sum()
    return int(counts)


def calculate_effint(counts, events):
    # Calculate efficiency and intensity
    effint = counts / events

    # Use binomial instead of poisson approximation
    effint_unc = np.sqrt(counts * (1 - counts/events)) / events

    return float(effint), float(effint_unc)


def calculate_B(counts, events, pseudo_time, measurement_time_hours):
    # Convert to seconds
    measurement_time = measurement_time_hours * 3600
    # Background count rate
    b = counts / pseudo_time
    b_unc = np.sqrt(counts * (1 - counts/events)) / pseudo_time
    # Background counts during measurement time
    B = b * measurement_time
    B_unc = b_unc * measurement_time
    return float(B), float(B_unc)


def calculate_B_filter(counts, events, Bq_, measurement_time_hours):
    # Activity of radionuclide in filter
    Bq, Bq_unc = Bq_

    # Convert to seconds
    measurement_time = measurement_time_hours * 3600
    # Calculate efficiency and intensity
    effint = counts / events
    # Use binomial instead of poisson approximation
    effint_unc = np.sqrt(counts * (1 - counts/events)) / events
    # Background counts
    B_filter = effint * Bq * measurement_time
    B_filter_unc = np.sqrt(np.power(Bq*measurement_time*effint_unc, 2) + np.power(effint*measurement_time*Bq_unc, 2))

    return B_filter, B_filter_unc


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


# TODO implement combined mda uncertainty calculation
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

