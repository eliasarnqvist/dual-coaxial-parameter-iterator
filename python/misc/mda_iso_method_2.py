import numpy as np
import matplotlib.pyplot as plt
# from scipy.stats import binom, poisson, beta, gamma
from scipy.optimize import brentq

np.random.seed(42)

# General
t_M = 3600 * 24 * 1
# Radionuclide
c_radionuclide = 1e4
N_radionuclide = 1e7
# Background SURE
c_SURE = 2e0
N_SURE = 1e9
t_SURE = 2e5
# Background filter
A = 0.200
A_unc = 0.020
c_filter = 1e3
N_filter = 1e7
# Cosmic
c_cosmic = 1e2
N_cosmic = 1e7
t_cosmic = 1e1
# For characteristic limits
alpha = 0.05
beta = 0.05

# Generate set of values for the background mean
N_rand = int(1e4)
# SURE
SURE_counts_mean_samples = np.random.gamma(c_SURE, scale=1, size=N_rand)
# SURE_counts_samples = np.random.poisson(lam=SURE_counts_mean_samples, size=N_rand)
b_SURE_samples = SURE_counts_mean_samples / t_SURE
# Filter
filter_effint_counts_mean_samples = np.random.gamma(c_filter, scale=1, size=N_rand)
# filter_effint_counts_samples = np.random.poisson(lam=filter_effint_counts_mean_samples, size=N_rand)
filter_effint_samples = filter_effint_counts_mean_samples / N_filter
filter_A_samples = np.random.normal(A, A_unc, size=N_rand)
b_filter_samples = filter_A_samples * filter_effint_samples

# background together
n_B_mean_samples = t_M * (b_SURE_samples + b_filter_samples)
n_B_samples = np.random.poisson(lam=n_B_mean_samples, size=N_rand)

y_means = np.array([])
for n_B_mean in n_B_mean_samples:
    N_rand_2 = int(1e2)

    radionuclide_effint_counts_mean_samples = np.random.gamma(c_radionuclide, scale=1, size=N_rand_2)
    radionuclide_effint_counts_samples = np.random.poisson(lam=radionuclide_effint_counts_mean_samples, size=N_rand_2)
    effint_samples = radionuclide_effint_counts_samples / N_radionuclide

    n_B_samples = np.random.poisson(lam=n_B_mean, size=N_rand_2)
    n_G_samples = np.random.poisson(lam=n_B_mean, size=N_rand_2)

    y_samples = 1 / (effint_samples * t_M) * (n_G_samples - n_B_samples)

    y_star = np.quantile(y_samples, 1 - alpha)

    def objective_function(n_G_value):
        n_G_value_samples = np.random.poisson(n_G_value, N_rand_2)
        y_hash_samples = 1 / (effint_samples * t_M) * (n_G_value_samples - n_B_samples)
        fraction_under_y_star = np.mean(y_hash_samples < y_star)
        return fraction_under_y_star - beta
    n_G_hash = brentq(objective_function, 1*n_B_mean, 10*n_B_mean)
    n_G_hash_samples = np.random.poisson(n_G_hash, N_rand_2)
    y_hash_samples = 1 / (effint_samples * t_M) * (n_G_hash_samples - n_B_samples)
    y_hash = np.mean(y_hash_samples)

    y_mean = np.mean(y_samples)
    y_means = np.append(y_means, y_hash)





# Plot
inch_to_mm = 25.4

fig, ax = plt.subplots(1, 1, figsize=(88/inch_to_mm, 60/inch_to_mm))

histo, ex = np.histogram(y_means, bins=100)
ax.step(ex[:-1], histo, where="post")
# histo, ex = np.histogram(n_B_samples2, bins=100)
# ax.step(ex[:-1], histo, where="post")

# ax.axvline(x=np.mean(SURE_counts_mean_samples), ymax=0.9, color='blue', ls='--')
# print(np.mean(SURE_counts_mean_samples))

ax.set_xlabel("y (cps)", fontsize=8)
ax.set_ylabel("pdf", fontsize=8)
plt.tight_layout(pad = 0.2)
fig.subplots_adjust(hspace=0, wspace=0)

plt.show()