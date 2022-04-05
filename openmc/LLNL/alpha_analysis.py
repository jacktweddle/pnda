import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

data = pd.read_excel('results_summary.xlsx', sheet_name='2.5 MeV, Isotropic - w TSL', names=('radius', 'buckling', 'published_alpha', 'exp_alpha', 'alpha', 'alpha_unc'), usecols=(0, 1, 3, 4, 5, 6), skiprows=3, nrows=14)

print(data)

v = 2.2E5  # Thermal neutron velocity in cm/s


def fit(x, D0, sigma_a, C):
    return v*sigma_a + D0*x - C*x**2


x_vals = np.linspace(0, 1, 1000)
popt, pcov = curve_fit(fit, data['buckling'], data['alpha']*1000, sigma=data['alpha_unc']*1000) #, p0=(26000, 0.033, 2500))
plt.plot(data['buckling'], data['published_alpha']*1000)
plt.errorbar(data['buckling'], data['alpha']*1000, capsize=6, linestyle='none', marker='o', markersize=3, color='black', ecolor='black')
#plt.plot(x_vals, fit(x_vals, popt[0], popt[1], popt[2]))
plt.title('α as a function of buckling - 20 x 1E8 particles')
plt.xlabel('Buckling [cm$^{-2}$]')
plt.ylabel('α [s$^{-1}$]')
plt.legend(('Published', 'OpenMC'))
plt.tick_params(direction='in',top='on',bottom='on',left='on',right='on')

print('D0 = {:.4f} +- {:.4f}'.format(popt[0],np.sqrt(pcov[0,0])))
print('sigma_a = {:.7f} +- {:.7f}'.format(popt[1],np.sqrt(pcov[1,1])))
print('C = {:.3f} +- {:.3f}'.format(popt[2],np.sqrt(pcov[2,2])))
