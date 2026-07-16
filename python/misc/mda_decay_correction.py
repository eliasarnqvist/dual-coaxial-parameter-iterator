import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq
import matplotlib as mpl
mpl.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman"]
})

def solve_u(m):
    # Solve for the value of u=lambda*t for a given m=b/lambda
    # Here, lambda is the decay constant and b is the background count rate

    f = lambda u: u - np.log(1 + (2.71/2.325)*np.sqrt(u/m) + 2*u)

    lo = 1e-12
    hi = 10.0

    # Ensure there are opposite signs
    if not f(lo) < 0:
        raise ValueError
    while f(hi) < 0:
        hi *= 2

    return brentq(f, lo, hi)

# m = np.logspace(np.log10(4e-2), np.log10(4e7), 100)
m = np.logspace(np.log10(1e-4), np.log10(1e6), 100)
u = np.array([solve_u(mi) for mi in m])

# Note that: m=t12*b/ln(2)
t12_b = m * np.log(2)
# Note that: u=n*ln(2)
n = u / np.log(2)

inch_to_mm = 25.4

fig, ax = plt.subplots(1, 1, figsize=(88/inch_to_mm, 60/inch_to_mm))
ax.plot(t12_b, n)
ax.set_xscale("log")
# ax.set_yscale("log")
ax.set_ylabel("$n$")
ax.set_xlabel("$b \, t_{1/2}$")

plt.tight_layout(pad = 0.2)
fig.subplots_adjust(hspace=0, wspace=0)
# save_name = "measured_spectra_comparison"
# plt.savefig(f'figures/{save_name}.jpg', dpi=600)
# plt.savefig(f'figures/{save_name}.pdf')


plt.show()

# https://www.geogebra.org/calculator/ttsqwn9b