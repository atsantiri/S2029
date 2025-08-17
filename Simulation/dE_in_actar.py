

import matplotlib.pyplot as plt
import numpy as np

beams = [['17F','3.4MeVu'], ['4He','1MeV'],['4He','5MeV'],['1H','800keV']]
c = ['b','g','magenta','orange']
lnst=['solid','dashed']

for i, b in enumerate(beams):
   x,dE,dE_rec = np.loadtxt(f'SRIM/{b[0]}_{b[1]}_ACTARgas_700mbar_IONIZ.txt',skiprows=28).T
   plt.plot(x*1e-7,dE*1e4, label=f'{b[0]}@{b[1]}', color=c[i])
plt.legend()  
plt.xlabel('x (mm)')
plt.ylabel('dE/dx (keV/mm)')    
plt.show()
