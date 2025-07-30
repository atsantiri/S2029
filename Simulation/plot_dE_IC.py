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
        axs[0].plot(ene, y/max(y), label=f'{b}@{p}mbar',color=c[i],linestyle=lnst[j])

axs[0].set_xlabel('Energy (MeV/u)')
axs[0].set_ylabel('Counts')
# axs[0].set_ylim(0,4.7e15)
axs[0].set_title('For 4.5 MeV/u in IC')

trim_energy=[]
with open('SRIM/17F_IC_70mbar_TRIM.txt') as f:
    for line in f:
        if line.startswith('T'):
            parts = line.split()
            trim_energy.append(float(parts[3])*1e-6/17) # TRIM output in eV, convert to MeV/u

axs[0].hist(trim_energy,bins=70,color='skyblue',alpha=0.5,edgecolor='black',label='TRIM 17F@70mbar', weights=np.ones_like(trim_energy) /250)

trim_energy_17o=[]
with open('SRIM/17O_IC_70mbar_TRIM.txt') as f:
    for line in f:
        if line.startswith('T'):
            parts = line.split()
            trim_energy_17o.append(float(parts[3])*1e-6/17) # TRIM output in eV, convert to MeV/u

axs[0].hist(trim_energy_17o,bins=50,color='magenta',alpha=0.5,edgecolor='black',label='TRIM 17O@70mbar', weights=np.ones_like(trim_energy_17o) /250)


axs[0].legend(loc='upper left')

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