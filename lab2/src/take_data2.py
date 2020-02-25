import pickle
import glob, os
# I need to switch to numpy fft
import time


def d2():
    """
    Abbreviation function: dual get hard-coded for week 1 data collection.
    Acquire two streams of data through pico sampler channels A and B.
    50 mV range, divisor=1
    """
    pack = ugradio.pico.capture_data('50mV', divisor=1, nsamples=16000, dual_mode=True)
    pack.shape = (2, -1, 16000)
    return pack

def d200_broken():
    data_chunk = []
    for i in range(100):
        data_chunk.append(ugradio.pico.capture_data('50mV', divisor=1, nblocks=100))
        print(str(i + 1) + '% complete')
    return data_chunk

# recommended: loops=100
def d200_new(loops):
    data_chunk = []
    for i in range(loops):
        data_chunk.append(ugradio.pico.capture_data('50mV', divisor=1, dual_mode=True, nsamples=42000,nblocks=100))
        print(str(i + 1) + '% complete')
    return data_chunk

def d100():
    time.sleep(120)
    data_chunk = []
    for i in range(100):
        data_chunk.append(ugradio.pico.capture_data('50mV', divisor=6, dual_mode=True))
        print(str(i + 1) + '% complete')
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


def unpickle_folder():
    re_block, im_block = [], []
    for file in glob.glob("./*"):
        data = pickle.load(open(file, "rb"))
        re_block.append(data['real'])
        im_block.append(data['imaginary'])
    return re_block, im_block