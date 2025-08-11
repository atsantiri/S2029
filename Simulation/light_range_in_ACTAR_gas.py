import pandas as pd
import matplotlib.pyplot as plt


types = [['alpha', '1MeV'],
         ['alpha', '5MeV'],
         ['p','800keV']]

for particle in types:
    df = pd.read_csv('./SRIM/RANGE_3D_{}_{}_ACTAR_700mbar.txt'.format(particle[0], particle[1]),
                     skiprows=20, delim_whitespace=True, header=None, encoding='latin1')
    range = df.iloc[:, 1]*1e-7 # Angstrom to mm
    plt.hist(range, bins=100, alpha=0.5, label='{} {}'.format(particle[1], particle[0]))
    print(f'Mean range for {particle[1]} {particle[0]}: {range.mean():.2f} mm')

plt.legend()
plt.show()