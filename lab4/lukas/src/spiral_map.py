from astropy.io import fits
import numpy as np

HI_rest = 1420.405751786e6 # MHz
c = 3e10 # cm / s

def frame():
    fig = plt.figure(figsize = (6, 3))
    plt.subplots_adjust(left=.15, bottom=.15, right=.95, top=.9)
    ax = fig.add_subplot(111)
    ax.tick_params(axis="x", labelsize=12)
    ax.tick_params(axis="y", labelsize=12)
    return fig, ax

def full_cal_plot(label, lon):
    fig, ax = frame()
    
    y, dopc = full_calibration(label, lon)
    x = np.linspace(1415e6, 1425e6, 8192) / 1e9

    plt.xlabel('RF Frequency [GHz]', fontsize=12)
    plt.ylabel('$T_{sys} + T_{ant, HI}$ [K]', fontsize=12)

    plt.plot(x, y)
    plt.vlines(HI_rest / 1e9, np.nanmin(y), np.nanmax(y))
    plt.show()

def doppler_plot(label, lon):
    fig, ax = frame()
    
    y, dopc = full_calibration(label, lon)
    frq = np.linspace(1415e6, 1425e6, 8192) / 1e9

    # Convert Doppler correction to km / s:
    dopc /= 1000

    x = [dopc + (1 - f / (HI_rest / 1e9)) * c / 1e5 for f in frq]

    plt.xlabel('Doppler Velocity [km / s]', fontsize=12)
    plt.ylabel('$T_{sys} + T_{ant, HI}$ [K]', fontsize=12)

    plt.plot(x, y)

    plt.vlines(dopc, np.nanmin(y), np.nanmax(y))
    plt.vlines(x[np.nanargmax(y)], np.nanmin(y), np.nanmax(y), color='orange')

    print(x[np.nanargmax(y)])

    plt.show()

# doppler_plot('cycle3', 120)

def calibration_fan(label, start_lon, stop_lon):
    return np.array([
        full_calibration(label, i) for i in range(start_lon, stop_lon + 1)
    ])

def full_calibration(label, lon):
    s_on_q, s_on_n, s_off_q, s_off_n, dopc = spectral_fan(label, lon)

    # T_noise is 270K for auto1, 80K for auto0
    gain_on = gain(s_on_n, s_on_q, 80)
    gain_off = gain(s_off_n, s_off_q, 80)
    gain_avg = .5 * (gain_on + gain_off)

    return gain_avg * s_on_q / s_off_q, dopc

    #we want to average each .fits file (ten spectra) into a single spectrum

def spectral_fan(label, angle, polarization = 0):
    '''
    Automatically expand a single reference into four file
    names for calibration:
    ex ffan('cycle_auto', 250)
    '''
    n = '_noisy.fits'
    nn = '_quiet.fits'
    lo1 = '_634MHz'
    lo2 = '_635MHz'
    pre = label + '_' + str(float(angle)) + '_degrees'

    s_on_quiet = read_average(pre + lo1 + nn, polarization)
    s_on_noisy = read_average(pre + lo1 + n, polarization)
    s_off_quiet = read_average(pre + lo2 + nn, polarization)
    s_off_noisy = read_average(pre + lo2 + n, polarization)

    # We arbitrarily take a correction at a place
    # which we hope is roughly in the middle of the observation.
    dopc = doppler_correction(pre + lo1 + nn)
    
    return s_on_quiet, s_on_noisy, s_off_quiet, s_off_noisy, dopc

def read_average(path, polarization, N=10):
    f = fits.open(path)

    assert polarization == 1 or polarization == 0
    key = 'auto0_real'
    if polarization == 1:
        key = 'auto1_real'

    spectra = np.array([f[i + 1].data[key] for i in range(N)])
    
    return np.average(spectra, axis=0)

def doppler_correction(path):
    '''
    Returns doppler correction based on the right ascension, declination,
    and Julian date extracted from the header on the fits file located
    at @path.
    We return a pure scalar which corresponds to the correction in m / s.
    '''
    f = fits.open(path)
    ra = f[0].header['RA']
    dec = f[0].header['DEC']
    jd = f[0].header['JD']
    return ugradio.doppler.get_projected_velocity(
        ra, dec, jd, ugradio.leo.lat, ugradio.leo.lon
    ).to_value()
 
'''
Current problems:
    there are three separate clouds. Which one represents the true Doppler
    shift? Why not all of them? The problem is that we have an
    equation for V_Dopp as a function of R.
'''
