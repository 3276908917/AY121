from astropy.io import fits
import numpy as np

HI_regular = 1420.405751786e6 # MHz
c = 3e10 # cm / s

def calibration_fan(label, start_lon, stop_lon):
    return np.array([
        full_calibration(label, i) for i in range(start_lon, stop_lon + 1)
    ])

def full_calibration(label, lon):
    s_on_q, s_on_n, s_off_q, s_off_n = spectral_fan(label, lon)
    
    gain_on = gain(s_on_n, s_on_q)
    gain_off = gain(s_off_n, s_off_q)

    return gain_on / gain_off * s_on_q / s_off_q

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
    
    return s_on_quiet, s_on_noisy, s_off_quiet, s_off_noisy

def read_average(path, polarization, N=10):
    f = fits.open(path)

    assert polarization == 1 or polarization == 0
    key = 'auto0_real'
    if polarization == 1:
        key = 'auto1_real'

    spectra = np.array([f[i + 1].data[key] for i in range(N)])
    
    return np.average(spectra, axis=0)

# draft function for collecting points to map
blue_shift():
    # for each
    f = fits.open(path)
    ra = f[0].header['RA']
    dec = f[0].header['DEC']
    jd = f[0].header['JD']
    correction = ugradio.doppler.get_projected_velocity(
        ra, dec, jd, ugradio.leo.lat, ugradio.leo.lon
    )

    # do we want to use a maximizer function to locate our HI line?
    
    '''
    Current problems:
        there are three separate clouds. Which one represents the true Doppler
        shift? Why not all of them? The problem is that we have an
        equation for V_Dopp as a function of R.
    '''
