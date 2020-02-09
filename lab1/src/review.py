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
        times.append(i / v_s * 10 ** 6)

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

def plot(x, y, xnote, ynote, im=False):
    """
    General but highly volatile plotting function
    (because I do not like writing functions in the shell).    
    Mostly ad-hoc and designed for rapid plotting 
    of similar categories of data.
    """
    plt.figure(figsize=(3,3))
    plt.subplots_adjust(left=.2, bottom=.15, right=.95, top=.9)
    if im:
        plt.plot(x, y[0], label='real component')
        plt.plot(x, y[1], label='imaginary component')
        plt.legend()
    else:    
        plt.plot(x, y)   
    
    plt.xlabel(xnote, fontsize=12)
    plt.ylabel(ynote, fontsize=12)
	
    plt.show()

def plot_sample(arr, norm, i=0, f=None):
    """
    Produce and display a plot
    where y-values are determined by @arr
    and where x-valies are determined by the current array associated with 'times'

    @i and @f are optional parameters to specify a contiguous slice of data
    """
    if f is None:
        f = len(arr)
    plt.figure(figsize=(3,3))
    plt.subplots_adjust(left=.2, bottom=.15, right=.95, top=.9)    
    plt.plot(times[i:f], normalize(arr[i:f], norm))	    
    plt.xlabel(r'Time ($\mu$s)', fontsize=12)
    plt.ylabel('Voltage (mV)', fontsize=12)
	
    plt.show()

    # I am assuming that we are continuing to use sample frequency 6.25e6 Hz
def raw_power(sample, rate=6.25e6):
    """
    Unlabeled power spectrum plot (with arbitrary normalization)
    of @sample at sampling rate @rate
    and return Fourier transform as well as power array

    This function is intended only for on-the-fly eyeballing.
    For reporting, use the analysis.py counterpart power_plot
    """
    f = ugradio.dft.dft(sample, vsamp=rate)
    P = np.abs(f[1]) ** 2
    plt.plot(f[0], P)
    plt.show()
    return f, P
