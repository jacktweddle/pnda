import openmc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

sp = openmc.StatePoint('statepoint.20.h5')

tally = sp.get_tally(name='detector1')
df = tally.get_pandas_dataframe()

pd.options.display.float_format = '{:,.3f}'.format
print(df)

filtered_data = df[(df['time low [s]'] >= 1E-3) & (df['time low [s]'] <= 1.5E-3)]
filtered_time = np.array(filtered_data['time low [s]']*1000)
filtered_flux = np.array(filtered_data['mean'])

print(filtered_time)

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

filtered_time = smooth(filtered_time)
filtered_flux = smooth(filtered_flux)

low_unc = np.inf
low_rel_err = np.inf

data = pd.DataFrame(columns=['Start time [ms]', 'End time [ms]', 'Alpha'])


# Function to find the optimal fit by searching for fits that fit different 'optimal' criteria
for i in range(25, len(filtered_time)-75, len(filtered_time)-75):
    for j in range(i+50, len(filtered_time), 10):
        time_interval = filtered_time[25:j]
        flux_interval = filtered_flux[25:j]
        popt, pcov = curve_fit(fit, time_interval, flux_interval, maxfev=500000)
        timestep_data = {'Start time [ms]': time_interval[0], 'End time [ms]': time_interval[-1], 'Alpha': popt[1]}
        data = data.append(timestep_data, ignore_index=True)
        # Find alpha value with lowest uncertainty
        if np.sqrt(pcov[1, 1]) < low_unc and popt[1] > 5:
            low_unc = np.sqrt(pcov[1, 1])
            low_unc_phi_0 = popt[0]
            low_unc_alpha = popt[1]
            low_unc_rr = popt[2]
            low_unc_time = time_interval
            low_unc_flux = flux_interval
        # Find alpha value with lowest relative error
        if np.sqrt(pcov[1, 1])/popt[1] < low_rel_err and popt[1] > 5:
            low_rel_err = np.sqrt(pcov[1, 1])/popt[1]
            low_rel_err_unc = np.sqrt(pcov[1, 1])
            low_rel_err_phi_0 = popt[0]
            low_rel_err_alpha = popt[1]
            low_rel_err_rr = popt[2]
            low_rel_err_time = time_interval
            low_rel_err_flux = flux_interval

plt.scatter(data['Start time [ms]'], data['End time [ms]'], c=data['Alpha'], vmin=7, vmax=10) #, vmax=15) #, norm=colors.LogNorm())
# plt.xlim(0, 1.05E-3)
plt.xlabel('Start time [ms]')
plt.ylabel('End time [ms]')
plt.colorbar(label='α [ms$^{-1}$]')
plt.show()

print('Lowest uncertainty')
print('Alpha =', low_unc_alpha, '±', low_unc, '\n')
print('Start time = ', low_unc_time[0])
print('End time = ', low_unc_time[-1])
print('######################################### \n')
print('Lowest relative error')
print('Alpha =', low_rel_err_alpha, '±', low_rel_err_unc, '\n')
print('Start time = ', low_rel_err_time[0])
print('End time = ', low_rel_err_time[-1])


plt.plot(df['time low [s]']*1000, np.abs(df['mean'])) #, s=0.1)
# plt.plot(filtered_time*1000, filtered_flux)
plt.xlim(0.9, 2)
plt.yscale('log')
#plt.ylim(1e-7, 1e-4)
plt.xlabel('Time [ms]')
plt.ylabel('Neutron current [particle/source]')
plt.title('Neutron current at detector - 0.66 ms pulse, water')
# plt.legend(labels=('Fit', 'Data', 'Filtered data'))
plt.show()


t_vals = np.linspace(1, 1.5, 1000)
popt, pcov = curve_fit(fit, filtered_time, filtered_flux) #, p0=(0.015, 35, 2.25E-5, 0.09, 0), bounds=((0, 30, 0, 0, 0), (0.1, 40, 1E-4, 0.2, 0.1)))
#plt.plot(t_vals, fit(t_vals, popt[0], popt[1], popt[2], popt[3], popt[4]), linestyle='-', color='black', linewidth=4)

plt.scatter(df['time low [s]']*1000, df['mean'], s=0.1)
plt.plot(t_vals, low_unc_phi_0*np.exp(-low_unc_alpha*t_vals))
#plt.plot(filtered_time*1000, filtered_flux)
plt.xlim(1, 1.2)
plt.ylim(0, 1.5e-6)
plt.xlabel('Time [ms]')
plt.ylabel('Neutron current [particle/source]')
plt.title('Neutron current at detector - 1ms pulse')
#plt.legend(labels=('Fit', 'Data', 'Smoothed data'))

'''
print('phi_0 = {:.7f} +- {:.7f}'.format(popt[0], np.sqrt(pcov[0, 0])))
print('alpha = {:.3f} +- {:.3f}'.format(popt[1], np.sqrt(pcov[1, 1])))
print('Room return = {:.7f} +- {:.7f}'.format(popt[2], np.sqrt(pcov[2, 2])))

print('C = {:.7f} +- {:.7f}'.format(popt[3], np.sqrt(pcov[3, 3])))
print('Room return = {:.8f} +- {:.8f}'.format(popt[4], np.sqrt(pcov[4, 4])))

'''
