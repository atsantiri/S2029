import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import AutoMinorLocator


# mpl.rcParams.update({
#     "font.family": "STIXGeneral",
#     "mathtext.fontset": "stix",
#     "axes.titlesize":17,
#     "axes.labelsize":20,
#     "font.size":20
# })

mpl.rcParams.update({
    "font.family": "STIXGeneral",
    "mathtext.fontset": "stix",
    "font.size": 20,
    "axes.labelsize": 20,
    "axes.titlesize": 17,
    "xtick.labelsize": 15,
    "ytick.labelsize": 15,
    "xtick.major.size": 6,
    "ytick.major.size": 6,
    "xtick.minor.size": 3,
    "ytick.minor.size": 3,
    "xtick.direction": "in",
    "ytick.direction": "in",
})

fname = '150-160deg'

x,deg,y,yerr = np.loadtxt(f'Outputs/doubleXS_{fname}.dat').T
ecn=x+3.92
plt.errorbar(ecn, y, yerr=yerr, fmt='.', capsize=3,color='black')
plt.xlabel(r'$E_{^{18}\mathrm{{Ne}}}$ [MeV]')
plt.ylabel(r'$\frac{d\sigma}{d\Omega}$ [$\frac{\mathrm{mb}}{\mathrm{sr}}$]')
plt.text(4.5,75,fr'$\theta_{{\mathrm{{CM}}}}$ = {deg[0]:.0f}$^o$ $\pm$ 5$^o$')

ax = plt.gca()
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())

ax.tick_params(which='both', direction='in')

plt.tight_layout()
plt.savefig(f'Outputs/xs_{fname}.png')
plt.show()
