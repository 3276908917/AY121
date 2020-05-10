# Functions to write:
    # 1. calculate reference frame adjustments
    # 2. average over frequency channels
        # "It's okay to degrade the frequency resolution to, say, 1 or 2 kHz"

#complex_combine
def cc(real, imag):
    return real + 1j*imag

def freq_range(v_s, N, W=1):
    '''
    Return a N-length array
    Where frequencies range between plus or minus W * v_s / 2

    CAUTION: this returns frequencies in units of Hertz!
        divide by 10 ** 6 if you want MHz.
    '''
    return np.linspace(- W * v_s / 2, W * v_s / 2)
    
def samples_median(power_list):
    '''
    Reduces @power_list (an array of power spectra)
    to a single power spectrum.
    '''
    return np.median(power_list, axis=0)

def samples_mean(power_list):
    '''
    Reduces @power_list (an array of power spectra)
    to a single power spectrum.
    '''
    return np.mean(power_list, axis=0)

def gain(scal, scold, T_sys=300):
    '''
    Return gain based on
    @scal :  power spectrum for thermal contamination
    @scold : power spectrum for cold sky data

    This is a direct codification of the lab instructions on gain
    calculation; refer to the lab for detailed explication of parameters.
    '''
    return sum(scold) * Tsys / sum(scal - scold)

# Doppler velocity in megameters per second
def doppler(nu, nu_0, offset):
    '''
    @nu_0: accepted value for HI line?
    @nu : current frequency in iteration over power spectrum

    This is a direct codification of the lab instructions on gain
    calculation; refer to the lab for explication of parameters.
    '''
    return -3e2 * (nu - nu_0) / nu_0 + offset

# deprecated
def power_lists(P, N):
    '''
    combines power level arrays, chiefly for use with the subsequent
    average and median methods.
    '''
    stack = [P[i][1] for i in range(1, N)]
    
