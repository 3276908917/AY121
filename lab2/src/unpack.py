import pickle
import glob, os

# Change these values as you run

# How many data points are there per axis (i.e. real + imaginary = 2)?
# If you have real and imaginary data: what is half the length of a block?
offset = 1600000

# How many data points are there per axis per sample?
# By default, the pico-sampler captures 16000 points per sample.
sample_size = 16000

# It's called bit because we are operating on a single block
# i.e. pico-sampler concatenates samples into one giant array.
def organize_bit(chunk):
    '''
    It is called bit because it is to be used on one block at a time.

    Return an array, equal in length to the number of blocks,
    of samples (arrays of complex voltages).

    For customization of the functionality, see global parameters
    in the header.
    '''
    c = []
    for i in range(0, offset, sample_size):
        real = chunk[i: i + sample_size]
        imag = chunk[offset + i: offset + i + sample_size]
        c.append(real + 1j*imag)
    return c

# e.g. off = load_saves('../on.npz')['raw_off']
def reduce_raw(case, reduction):
    '''
    Given a collection of blocks (pico-sampler concatenated samples),
    we turn each block into a group of arrays (in complex-conjugate
    format), then we perform the power spectrum on all arrays in that
    block, then we perform the reduction algorithm on each group
    of complex arrays, then we perform the reduction on all blocks.

    The two options for @reduction,
        samples_median and samples_mean,
    are described in data_appraisal.py

    Be careful about using this function in low-memory environments!
    The shell would either kill or be killed by the OS.
    '''
    compressor = []
    for a in case:
        c_chunk = organize_bit(a)
        P = power_barrage(c_chunk)
        compressor.append(reduction(P))
    return reduction(compressor)

# This helped me with Kyle's data.
def thermal():
    '''
    Assumes all files in current directory are data files,
        of the form: 1 float per line
    Copies the directory into an array, where each element
        is an array describing the floats contained within one file.
    '''
    block = []
    for file in glob.glob("./*"):
        samples = []
        reader = open(file, 'r')
        for line in reader:
            trim = line.strip()
            samples.append(float(trim))
        block.append(samples)
    return block

# This helped me with Max's data.
# It may help again in the future, in case I shift to pickle format.
def unpickle_folder():
    re_block, im_block = [], []
    for file in glob.glob("./*"):
        data = pickle.load(open(file, "rb"))
        re_block.append(data['real'])
        im_block.append(data['imaginary'])
    return re_block, im_block
