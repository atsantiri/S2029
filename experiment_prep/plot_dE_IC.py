# macro to plot the energy loss in the ionization chamber from LISE++ output files
# Looking at 20Mg and 20Na at 6MeV/u and 17F and 17O at 4.5MeV/u
# All beams pass through 11.5um mylar and 36mm CF4 at 70mbar

import matplotlib.pyplot as plt
import numpy as np

beams1 = ['17F', '17O']
beams2 = ['20Na', '20Mg']

fig, axs = plt.subplots(1,2, figsize=(10, 6))

for b in beams1:
    with open(f'LISEpp_outputs/{b}_70mbar_IC.txt') as f:
        clean_lines = (line.replace('-', '0') for line in f)
        ene, y = np.loadtxt(clean_lines, skiprows=4).T
    axs[0].plot(ene, y, label=b)

axs[0].legend()
axs[0].set_ylabel('Counts')
axs[0].set_title('For 4.5 MeV/u in 70mbar IC')

for b in beams2:
    with open(f'LISEpp_outputs/{b}_70mbar_IC.txt') as f:
        clean_lines = (line.replace('-', '0') for line in f)
        ene, y = np.loadtxt(clean_lines, skiprows=4).T
    axs[1].plot(ene, y, label=b)

axs[1].legend()
axs[1].set_xlabel('Energy (MeV/u)')
axs[1].set_ylabel('Counts')
axs[1].set_title('For 6 MeV/u in 70 mbar IC')

plt.tight_layout()
plt.show()