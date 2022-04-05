import openmc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.colors as colors
from scipy.optimize import curve_fit

sp = openmc.StatePoint('statepoint.10.h5')

tally = sp.get_tally(name='detector')

df = tally.get_pandas_dataframe()

# df.to_excel('nile_mirrobor_box.xlsx')

pd.options.display.float_format = '{:,.3f}'.format
print(df)

filtered_data = df[(df['time low [s]'] > 1.05E-3) & (df['time low [s]'] < 3E-3)]
filtered_time = np.array(filtered_data['time low [s]'])
filtered_flux = np.array(filtered_data['mean'])

v = 2.2E5  # thermal neutron velocity in cm/s
D = 0.16  # Diffusion coefficient for thermalised light water medium
delta = 2*D  # Extrapolation distance in diffusion theory - http://www.gammaexplorer.com/wp-content/uploads/2014/03/Introduction-to-Nuclear-Engineering-Lamarsh-3rd-Edition.pdf P254
h = 17.9  # Moderator height in cm
r = 8.9  # Moderator radius in cm
buckling = (np.pi/(h+(2*delta)))**2 + (2.405/(r+delta))**2


"""def fit(t, phi_0, D_0, sig_a, C, R):
    return phi_0*np.exp(-((v*sig_a)+(buckling*D_0)-(buckling**2*C))*t) + R"""


def fit(t, phi_0, alpha, R):
    return phi_0*np.exp(-alpha*t) + R


def smooth(data):
    """5 point smoothing function.
    Recommended for use in low statistics in
    G.W. Phillips , Nucl. Instrum. Methods 153 (1978), 449
    Parameters
    ----------
    """
    smooth_spec = []
    smooth_spec.append(data[0])
    smooth_spec.append(data[1])

    spec_len = len(data)
    i = 2
    while i < spec_len - 2:
        val = (1.0 / 9.0) * (
            data[i - 2]
            + data[i + 2]
            + (2 * data[i + 1])
            + (2 * data[i - 1])
            + (3 * data[i])
        )
        smooth_spec.append(val)
        i = i + 1
    smooth_spec.append(data[i])
    smooth_spec.append(data[i + 1])

    return np.array(smooth_spec)


# To smooth filtered data
# filtered_time = smooth(filtered_time)
# filtered_flux = smooth(filtered_flux)

# To fit to unsmoothed data
#filtered_time = np.array(filtered_data['time low [s]'])
#filtered_flux = np.array(filtered_data['mean'])

plt.scatter(df['time low [s]'], df['energy low [eV]'], c=df['mean'], vmin=0) #, vmax=1e-7) #, norm=colors.LogNorm())
# plt.xlim(0, 1.05E-3)
plt.xlabel('Time [ms]')
plt.ylabel('Energy [eV]')
plt.colorbar(label='Neutron current [particle/source]')
plt.show()


'''
t_vals = np.linspace(1.05, 1.5, 1000)
plt.plot(df['time low [s]']*1000, df['mean'])
popt, pcov = curve_fit(fit, filtered_time*1000, filtered_flux) #, p0=(0, 8, 0), bounds=((0, 6, 0), (1, 11, 1)))
plt.plot(t_vals, fit(t_vals, popt[0], popt[1], popt[2]), linestyle='-', color='black', linewidth=3)
plt.plot(filtered_time*1000, filtered_flux)

plt.xlim(0.9, 1.5)
plt.xlabel('Time [ms]')
plt.ylabel('Neutron current [particle/source]')
plt.title('Neutron current at detector - 1 ms pulse, poly')
plt.legend(labels=('Fit', 'Data', 'Filtered data'))

print('phi_0 = {:.7f} +- {:.7f}'.format(popt[0], np.sqrt(pcov[0, 0])))
print('alpha = {:.3f} +- {:.3f}'.format(popt[1], np.sqrt(pcov[1, 1])))
print('Room return = {:.7f} +- {:.7f}'.format(popt[2], np.sqrt(pcov[2, 2])))
"""
print('C = {:.7f} +- {:.7f}'.format(popt[3], np.sqrt(pcov[3, 3])))
print('Room return = {:.8f} +- {:.8f}'.format(popt[4], np.sqrt(pcov[4, 4])))
"""
'''
