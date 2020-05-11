from astropy.io import fits
import numpy as np

HI_rest = 1420.405751786e6 # MHz
c = 3e10 # cm / s

V0 = 220 # km / s
R0 = 2.62282594e17 # km
kpc = 3.08567758e16 # km to kpc conversion

def frame():
    '''
    Set up generic infrastructure for an improved-looking plot.
    We return the figure and the axis on which we want to plot.
    '''
    fig = plt.figure(figsize = (6, 3))

    plt.subplots_adjust(left=.15, bottom=.15, right=.95, top=.9)
    ax = fig.add_subplot(111)
    
    ax.tick_params(axis="x", labelsize=12)
    ax.tick_params(axis="y", labelsize=12)

    return fig, ax

def arrow(doppler):
    '''
    Calculate a maximum-breadth velocity curve for a single
    galactic longitude and its associated Doppler velocity.
    @doppler is in the form (ell, Doppler_velocity)
    '''
    sin_ell = np.sin(np.radians(doppler[0]))
    #v_dopp = doppler[1]
    radial_space = np.linspace(R0 * abs(sin_ell), R0, 50)
    #velocity = lambda r: r / R0 * (v_dopp / sin_ell - V0)

    velocity = velocify(doppler)

    return radial_space, np.array([velocity(r) for r in radial_space])

def full_yoke(dopplers):
    '''
    For every doppler pair (ell, Doppler_velocity)
    in the list @dopplers,
    plot the associated arrow.
    We plot everything at once.
    '''
    fig, ax = frame()

    plt.xlabel('Distance from Galactic Center [kpc]', fontsize=12)
    plt.ylabel('Linear Velocity [km / s]', fontsize=12)

    for datum in dopplers:
        r, v = arrow(datum)
        plt.plot(r / kpc, v)

    plt.show()

def skip_yoke(dopplers):
    '''
    We take in a list of doppler pairs @dopplers,
    where each entry is a pair in the form (ell, Doppler_velocity),
    and we plot one arrow for every ten degrees of longitude.
    We plot everything at once, but because we are yoking
    fewer total arrows, we can maintain a legend entry
    describing the angle corresponding to each arrow.
    '''
    fig, ax = frame()

    plt.xlabel('Distance from Galactic Center [kpc]', fontsize=12)
    plt.ylabel('Linear Velocity [km / s]', fontsize=12)

    for i in range(0, len(dopplers), 10):
        r, v = arrow(dopplers[i])
        plt.plot(r / kpc, v, label=str(dopplers[i][0]) + ' degrees')

    #ax.legend(loc='lower left')
    ax.legend(bbox_to_anchor=(1, 1))
    
    plt.show()

def full_cal_plot(label, lon):
    '''
    Given the label of the data @label
    and the corresponding value of galactic longitude @lon
    Plot the fully calibrated and labeled spectrum as per
    equation \ref{eq:line_shape} in my lab 4 report WIP.

    You have to be sitting in the data directory for this to work.
    '''
    fig, ax = frame()
    
    y, dopc = full_calibration(label, lon)
    x = np.linspace(1415e6, 1425e6, 8192) / 1e9

    plt.xlabel('RF Frequency [GHz]', fontsize=12)
    plt.ylabel('$T_{sys} + T_{ant, HI}$ [K]', fontsize=12)

    plt.plot(x, y)
    plt.vlines(HI_rest / 1e9, np.nanmin(y), np.nanmax(y))
    plt.show()

def doppler_fan(label, start_lon, stop_lon, show=False):
    '''
    Create an array of the LSR Doppler-corrected fully
    calibrated spectra for every galactic longitude
    starting with @start_lon and ending with @stop_lon,
    where the files have @label in common.

    You have to be sitting in the data directory for this to work.
    '''
    doppler_collection = []
    for ell in range(start_lon, stop_lon + 1):
        if show:
            peak = full_doppler_plot(label, ell)
        else:
            x, y, dopc, peak = full_doppler(label, ell)
        doppler_collection.append((ell, peak))
    return np.array(doppler_collection)

def full_doppler_plot(label, lon):
    '''
    Given the label of the data @label
    and the corresponding value of galactic longitude @lon
    plot the LSR Doppler-corrected fully calibrated and labeled spectrum.
    Also return the value for the Doppler velocity at the point
    which is deemed to be the maximum. 
    
    You have to be sitting in the data directory for this to work.
    '''
    fig, ax = frame()

    plt.xlabel('Doppler Velocity [km / s]', fontsize=12)
    plt.ylabel('$T_{sys} + T_{ant, HI}$ [K]', fontsize=12)

    x, y, dopc, peak_freq = full_doppler(label, lon)
    plt.plot(x, y)

    plt.vlines(dopc, np.nanmin(y), np.nanmax(y))
    plt.vlines(x[np.nanargmax(y)], np.nanmin(y), np.nanmax(y), color='orange')

    print('Peak Velocity [km / s]', peak_freq)
    plt.show()
    return peak_freq

def full_doppler(label, lon):
    '''
    Given a data label @label and a galactic longitude @lon,
    compute the fully calibrated spectrum and represent the
    x-axis in Doppler velocities corrected for the LSR frame.

    You have to be sitting in the data directory for this to work.
    '''
    y, dopc = full_calibration(label, lon)
    frq = np.linspace(1415e6, 1425e6, 8192) / 1e9
    
    dopc /= 1000
    x = [dopc + (1 - f / (HI_rest / 1e9)) * c / 1e5 for f in frq]
    
    peak_freq = x[np.nanargmax(y)]
    return x, y, dopc, peak_freq

def calibration_fan(label, start_lon, stop_lon):
    '''
    Create an array of the fully calibrated spectra for every galactic longitude
    starting with @start_lon and ending with @stop_lon,
    where the files have @label in common.

    You have to be sitting in the data directory for this to work.
    '''
    return np.array([
        full_calibration(label, i) for i in range(start_lon, stop_lon + 1)
    ])

def full_calibration(label, lon):
    '''
    Given the label of the data @label
    and the corresponding value of galactic longitude @lon
    return the fully calibrated and labeled spectrum as per
    equation \ref{eq:line_shape} in my lab 4 report WIP.

    You have to be sitting in the data directory for this to work.
    '''
    s_on_q, s_on_n, s_off_q, s_off_n, dopc = spectral_fan(label, lon)

    # T_noise is 270K for auto1 (horizontal polarization, which we exclude),
    # 80K for auto0 (vertical polarization)
    gain_on = gain(s_on_n, s_on_q, 80)
    gain_off = gain(s_off_n, s_off_q, 80)
    gain_avg = .5 * (gain_on + gain_off)

    return gain_avg * s_on_q / s_off_q, dopc

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
    '''
    Return the average of the @N spectra contained
    in the .fits file located at @path
    and associated with the polarization @polarization.
    '''
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
