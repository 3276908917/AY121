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
    pack = ugradio.pico.capture_data('50mV', divisor=6, nsamples=16000, dual_mode=True)
    pack.shape = (2, -1, 16000)
    return pack

# recommended: loops=100
def d200(loops):
    data_chunk = []
    for i in range(loops):
        data_chunk.append(ugradio.pico.capture_data('50mV', divisor=6, dual_mode=True, nblocks=100))
        print(str(i + 1) + '% complete')
    return data_chunk

def complex_combine(real, imag=None):
    if imag is None:
        r = real[100:16000]
        i = real[16100:32000]
        return r + 1j*i
    return real + 1j*imag

def data_to_comp(glob, sample_size=16000):
    complex_combo = []    
    #samples_per_block = len(glob[0][0])
    offset = len(glob[0]) // 2
 #   num_pairs = offset / sample_size

    for c in glob:
        for i in range(0, offset, sample_size):
            real = c[i: i + sample_size]
            imag = c[offset + i: offset + i + sample_size]
            complex_combo.append(real + 1j*imag)

    return complex_combo

def complex_bblock(glob):
    c = []
    offset = len(glob[0]) // 2
    for a in glob:
        c.append(a[0:offset] + 1j*a[offset:])
    return c

# 100 blocks of 320000
