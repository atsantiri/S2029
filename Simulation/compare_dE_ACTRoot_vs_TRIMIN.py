# macro to plot the energy loss in the ionization chamber from LISE++ output files
# Looking at 20Mg and 20Na at 6MeV/u and 17F and 17O at 4.5MeV/u
# All beams pass through 11.5um mylar and 36mm CF4 at 70mbar

import matplotlib.pyplot as plt
import numpy as np

files = [['17F after IC TRIM','17F_IC_70mbar_TRIM.txt'],
         ['17O after IC TRIM','17O_IC_70mbar_TRIM.txt'],
         ['17F after CFA TRIM','17F_IC70mbar_CFA6mbar_TRIMIN.txt'],
         ['17O after CFA TRIM','17O_IC70mbar_CFA6mbar_TRIMIN.txt'],
         ]

c = ['b','g','magenta','orange']
lnst = ['solid', 'dashed']
fig, axs = plt.subplots(1, 1, figsize=(8, 6))

axs.set_xlabel('Energy (MeV/u)')
axs.set_ylabel('Counts')

# Plot TRIM histograms (solid lines)
for i, file in enumerate(files):
    energy = []
    with open(f'SRIM/{file[1]}') as f:
        for line in f:
            if line.startswith('T'):
                parts = line.split()
                energy.append(float(parts[3]) * 1e-6 / 17)  # TRIM output in eV, convert to MeV/u

    axs.hist(energy,bins=70,alpha=1.0,edgecolor=c[i],facecolor='none',label=file[0],density=True,linestyle=lnst[0],linewidth=1.5,histtype='step')

# Load ACTRoot data
actroot_17FIC, actroot_17OIC, actroot_17FCFA, actroot_17OCFA = np.loadtxt('dE_out.txt', skiprows=1).T
actroot_data = [
    (actroot_17FIC, '17F after IC ACTRoot', c[0]),
    (actroot_17OIC, '17O after IC ACTRoot', c[1]),
    (actroot_17FCFA, '17F after CFA ACTRoot', c[2]),
    (actroot_17OCFA, '17O after CFA ACTRoot', c[3]),
]

# Plot ACTRoot histograms (dashed lines)
for data, label, color in actroot_data:
    axs.hist(data,bins=50,alpha=1.0,edgecolor=color,facecolor='none',label=label,density=True,linestyle=lnst[1],linewidth=1.5,histtype='step')

plt.legend()

plt.tight_layout()
plt.savefig('compare_TRIM_ACTRoot.png')
plt.show()