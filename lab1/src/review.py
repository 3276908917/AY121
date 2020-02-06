import numpy as np
import matplotlib.pyplot as plt

# used to plot sampled waveforms; facilitates visual inspection
times = []

def create_time(v_s):
    """
    Rewrite the global variable times:
    length-16000 array of elapsed times,
    where each index gives a time according to
    @v_s sampling frequency.
    """
    global times
    for i in range(0, 16000):
        times.append(i / v_s)

# this default line poises the script for lab 1 week 1 data collection
create_time(6.25e6)

def load_saves(filename):
    """
    Return a dictionary containing the arrays saved to
    the .npz file at @filename
    """
    a = np.load(filename)
    return dict(zip(("{}".format(k) for k in a), (a[k] for k in a)))

def normalize(arr, actual_max):
    """
    Return a copy of @arr
    where each value has been scaled such that its max value
    corresponds to the given @actual_max
    """
    norm_c = actual_max / max(arr)
    return arr * norm_c

def plot(arr)
    """
    Produce and display a plot
    where y-values are determined by @arr
    and where x-valies are determined by the current array associated with 'times'
    """
    plt.plot(times[:len(arr)], arr)
    plt.show()
