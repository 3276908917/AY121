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

def doppler_fan(label, start_lon, stop_lon, show=False):
    doppler_collection = []
    for ell in range(start_lon, stop_lon + 1):
        if show:
            peak = full_doppler_plot(label, ell)
        else:
            x, y, peak = full_doppler_plot(label, ell)
        doppler_collection.append((ell, peak))
    return np.array(doppler_collection)

def full_doppler_plot(label, lon):
    fig, ax = frame()

    plt.xlabel('Doppler Velocity [km / s]', fontsize=12)
    plt.ylabel('$T_{sys} + T_{ant, HI}$ [K]', fontsize=12)

    x, y, peak_freq = full_doppler(label, lon)
    plt.plot(x, y)

    plt.vlines(dopc, np.nanmin(y), np.nanmax(y))
    plt.vlines(x[np.nanargmax(y)], np.nanmin(y), np.nanmax(y), color='orange')

    print('Peak frequency [Hz]', peak_freq)
    plt.show()
    return peak_freq

def full_doppler(label, lon):
    y, dopc = full_calibration(label, lon)
    frq = np.linspace(1415e6, 1425e6, 8192) / 1e9
    dopc /= 1000
    x = [dopc + (1 - f / (HI_rest / 1e9)) * c / 1e5 for f in frq]
    peak_freq = x[np.nanargmax(y)]
    return x, y, peak_freq

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
    
    
    We can get the doppler shift with respect to temperature by using
    1D interpolation scipy.interpolate.interp1d on the corrected dopper velocity 
    and have it fit the frequency function, which is also the same frequency function used with the temperatures.
    The frequency function for both velocity and temperature is based on the sampling rate of our telecope,
    using   t = 1/header["SAMPRATE"]
    
    It's basically graphing the velocity vs frequency and temperature vs frequency, to get at the end temperature vs velocity.
	
	We can use the Temperature vs doppler velocity function to calculate the HI density N_HI (function 3 in the lab handout)
    and Mass of HI M_Hi (function 7 in the lab handout) at any Galactic longtitude l, b=0. 
	
	To detect the spiral movements, We can make a 2d image using the Galactic longtitude as x axis, 
    Doppler Velocity as y axis and temperature as color. Grid[ Doppler velocity ][ Galactic longtitude ] = Temperature    
    and plot using plt.imshow(Grid)
	
	We can also create 2d image using the Galactic longtitude as x axis, density N of HI as y axis and Doppler velocity as input for color. 
    Grid[ Doppler velocity ][ Galactic longtitude ] = density N HI
    
    Then we can make a movie by ploting an image of 1 * 260 pixels where the y is Galactic Latitude = 0,
    y axis is Galactic longtitude from -10 to 250 degrees, 
    and the color shows how the density N_Hi is changing by moving the velocity from -200 km/s to 200 km/s
    
    
'''
