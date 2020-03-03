# vertical angle difficult to measure: plus or minus 10 degrees,
# which is already the angular resolution, probably of the big horn

import time

def d2():
    '''
    Abbreviation function: dual get hard-coded for week 1 data collection.
    Acquire two streams of data through pico sampler channels A and B.
    50 mV range, divisor=6
    '''
    pack = ugradio.pico.capture_data('50mV', divisor=6, nsamples=16000, dual_mode=True)
    pack.shape = (2, -1, 16000)
    return pack

# recommended: loops=100
def d200(loops):
    '''
    Abbreviation function: capture pico-sampler dual channel data,
    at 50 mV range, with divisor 6, where each increment of
    @loops corresponds to an additional order of 100 blocks of data.
    E.g. loops=100 -> 10k samples returned

    Progress bar is hard-coded for loops=100
    '''
    data_chunk = []
    for i in range(loops):
        data_chunk.append(ugradio.pico.capture_data('50mV', divisor=6, dual_mode=True, nblocks=100))
        print(str(i + 1) + '% complete')
    return data_chunk

def sleeper(seconds):
    time.sleep(seconds)    
    return ugradio.pico.capture_data('50mV', divisor=6, dual_mode=True, nblocks=100)

