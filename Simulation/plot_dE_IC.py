# macro to plot the energy loss in the ionization chamber from LISE++ output files
# Looking at 20Mg and 20Na at 6MeV/u and 17F and 17O at 4.5MeV/u
# All beams pass through 11.5um mylar and 36mm CF4 at 70mbar

import matplotlib.pyplot as plt
import numpy as np

beams1 = ['17F', '17O']
beams2 = ['20Na', '20Mg']
c = ['b','g','magenta']
lnst=['solid','dashed']
fig, axs = plt.subplots(1,2, figsize=(10, 6))
pressure = ['60','70','80']

for j,b in enumerate(beams1):
    for i,p in enumerate(pressure):
        with open(f'LISEpp_outputs/{b}_{p}mbar_IC.txt') as f:
            clean_lines = (line.replace('-', '0') for line in f)
            ene, y = np.loadtxt(clean_lines, skiprows=4).T
        axs[0].plot(ene, y, label=f'{b}@{p}mbar',color=c[i],linestyle=lnst[j])

axs[0].legend(loc='upper left')
axs[0].set_xlabel('Energy (MeV/u)')
axs[0].set_ylabel('Counts')
axs[0].set_ylim(0,4.7e15)
axs[0].set_title('For 4.5 MeV/u in IC')

for j,b in enumerate(beams2):
    for i,p in enumerate(pressure):
        with open(f'LISEpp_outputs/{b}_{p}mbar_IC.txt') as f:
            clean_lines = (line.replace('-', '0') for line in f)
            ene, y = np.loadtxt(clean_lines, skiprows=4).T
        axs[1].plot(ene, y, label=f'{b}@{p}mbar',color=c[i],linestyle=lnst[j])

axs[1].legend(loc='upper left')
axs[1].set_xlabel('Energy (MeV/u)')
axs[1].set_ylabel('Counts')
axs[1].set_ylim(0,2.5e15)

axs[1].set_title('For 6 MeV/u in IC')

plt.tight_layout()
plt.savefig('dE_IC_vary_pressure.png')
plt.show()