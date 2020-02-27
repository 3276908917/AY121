# Functions to write:
    # 1. calculate reference frame adjustments
    # 2. average over frequency channels
        # "It's okay to degrade the frequency resolution to, say, 1 or 2 kHz"

#complex_combine
def cc(real, imag):
    return real + 1j*imag

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

def gain(scal, scold):
    '''
    Return gain based on
    @scal :  power spectrum for thermal contamination
    @scold : power spectrum for cold sky data

    This is a direct codification of the lab instructions on gain
    calculation; refer to the lab for detailed explication of parameters.
    '''
    return sum(scold) * 300/sum(scal - scold)

# I guess the line frequency is the accepted value for the HI line?
def doppler(nu_0, nu):
    '''
    @nu_0
    @nu : current frequency in iteration over power spectrum

    This is a direct codification of the lab instructions on gain
    calculation; refer to the lab for explication of parameters.
    '''
    return 3e10 * (nu_0 - nu) / nu_0

# deprecated
def power_lists(P, N):
    '''
    combines power level arrays, chiefly for use with the subsequent
    average and median methods.
    '''
    stack = [P[i][1] for i in range(1, N)]
    
