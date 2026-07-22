import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import binom, beta
from scipy.optimize import brentq

# np.random.seed(42)

MDAs = []

random_tests = int(1e3)
random_seeds = [i for i in range(random_tests)]
for random_seed in random_seeds:
    np.random.seed(random_seed)

    # General
    t_M = 3600 * 24 * 1
    # Radionuclide
    c_radionuclide = 1e5
    N_radionuclide = 1e7
    p_radionuclide_estimator = c_radionuclide / N_radionuclide
    # Background SURE
    c_SURE = 1e1
    N_SURE = 1e9
    p_SURE_estimator = c_SURE / N_SURE
    t_SURE = 2e5
    # Background filter
    A = 0.300
    A_unc = 0.01
    c_filter = 1e2
    N_filter = 1e7
    p_filter_estimator = c_filter / N_filter
    # Cosmic
    c_cosmic = 1e2
    N_cosmic = 1e7

    # Generate random set of values
    N_rand = int(1e3)
    # For effint
    c_radionuclide_samples = np.random.binomial(N_radionuclide, p_radionuclide_estimator, size=N_rand)
    effint_radionuclide_samples = c_radionuclide_samples / N_radionuclide
    # For SURE
    c_SURE_samples = np.random.binomial(N_SURE, p_SURE_estimator, size=N_rand)
    b_SURE_samples = c_SURE_samples / t_SURE
    # For filter
    c_filter_samples = np.random.binomial(N_filter, p_filter_estimator, size=N_rand)
    effint_filter_samples = c_filter_samples / N_filter
    A_filter_samples = np.random.normal(A, A_unc, size=N_rand)
    b_filter_samples = A_filter_samples * effint_filter_samples
    # All together now
    b_samples = b_SURE_samples + b_filter_samples
    n_B_samples = t_M * b_samples
    n_B_mean = t_M * (c_SURE/t_SURE + A*p_filter_estimator)
    n_G_samples = np.random.poisson(n_B_mean, N_rand)

    y_samples = (1 / effint_radionuclide_samples) * (1 / t_M) * (n_G_samples - n_B_samples)

    # Determine the detection limit using the MC approach
    # False positive rate
    alpha = 0.05
    # False negative rate
    beta = 0.05
    # Decision threshold, critical limit, LC, y*
    y_star = np.quantile(y_samples, 1 - alpha)
    # Detection limit, LD, y#
    def objective_function(n_G_value):
        n_G_value_samples = np.random.poisson(n_G_value, N_rand)
        y_hash_samples = (1 / effint_radionuclide_samples) * (1 / t_M) * (n_G_value_samples - n_B_samples)
        fraction_under_y_star = np.mean(y_hash_samples < y_star)
        return fraction_under_y_star - beta
    n_G_hash = brentq(objective_function, 1*n_B_mean, 10*n_B_mean)
    n_G_hash_samples = np.random.poisson(n_G_hash, N_rand)
    y_hash_samples = (1 / effint_radionuclide_samples) * (1 / t_M) * (n_G_hash_samples - n_B_samples)
    y_hash = np.mean(y_hash_samples)

    # print(y_star, np.mean(y_samples > y_star))
    # print(y_hash, np.mean(y_hash_samples < y_star))

    MDAs.append(y_hash)

MDAs = np.array(MDAs)
MDA = np.mean(MDAs)

y_0_histo, ex_0 = np.histogram(y_samples, bins=100)
y_0_histo_norm = y_0_histo / N_rand

y_hash_histo, ex_hash = np.histogram(y_hash_samples, bins=100)
y_hash_histo_norm = y_hash_histo / N_rand

MDA_histo, ex_MDA = np.histogram(MDAs, bins=20)
MDA_mean = np.mean(MDAs)
MDA_std = np.std(MDAs)
print(MDA_mean, MDA_std)

# Plot
inch_to_mm = 25.4

fig, ax = plt.subplots(1, 1, figsize=(88/inch_to_mm, 60/inch_to_mm))
ax.step(ex_0[:-1], y_0_histo_norm, where="post")
ax.axvline(x=y_star, ymax=0.9, color='red', ls='--')
ax.step(ex_hash[:-1], y_hash_histo_norm, where="post")
ax.axvline(x=y_hash, ymax=0.9, color='blue', ls='--')
ax.set_xlabel("y (Bq)", fontsize=8)
ax.set_ylabel("pdf", fontsize=8)
plt.tight_layout(pad = 0.2)
fig.subplots_adjust(hspace=0, wspace=0)

fig, ax = plt.subplots(1, 1, figsize=(88/inch_to_mm, 60/inch_to_mm))
ax.step(ex_MDA[:-1], MDA_histo, where="post")
ax.set_xlabel("MDA (Bq)", fontsize=8)
ax.set_ylabel("distribution", fontsize=8)
plt.tight_layout(pad = 0.2)
fig.subplots_adjust(hspace=0, wspace=0)

plt.show()

