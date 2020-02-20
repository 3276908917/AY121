# I need to switch to numpy fft

def d2():
    """
    Abbreviation function: dual get hard-coded for week 1 data collection.
    Acquire two streams of data through pico sampler channels A and B.
    50 mV range, divisor=1
    """
    return ugradio.pico.capture_data('50mV', divisor=1, dual_mode=True)

def d200():
    data_chunk = []
    for i in range(10):
        for j in range(100):
            data_chunk.append(d2())
        print('10% increment')
    return data_chunk


def complex_combine(real, imag=None):
    if imag is None:
        r = real[100:16000]
        i = real[16100:32000]
        return r + 1j*i
    return real + 1j*imag

# combines power level arrays, chiefly for use with the subsequent
# average and median methods.
def power_lists(P, N):
    stack = [P[i][1] for i in range(1, N)]

def samples_median(power_list):
    return np.median(power_list, axis=0)

def samples_mean(power_list):
    return np.mean(power_list, axis=0)

# "It's okay to degrade the frequency resolution to, say, 1 or 2 kHz'
	# I should probably write a function to handle this.

# Big: calculate reference frame adjustments

def gain(scal, scold):
    return sum(scold) * 300/sum(scal - scold)

# I guess the line frequency is the accepted value for the HI line?
def doppler(nu_0, nu):
    return 3e10 * (nu_0 - nu) / nu_0
