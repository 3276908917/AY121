import numpy as np
import matplotlib.pyplot as plt

# used to plot sampled waveforms; facilitates visual inspection
times = []

def create_time(v_s, L=16000):
    """
    Rewrite the global variable times:
    length-16000 array of elapsed times,
    where each index gives a time according to
    @v_s sampling frequency.
    """
    global times
    for i in range(0, L):
        times.append(i / v_s * 10 ** 6)

# I use this for 5.4, but scale and N are not actually separable,
    # so this function may not be doing what I hope it to be doing...
def scale_freqs(v_s, N, scale):
    """
    Return a length @scale*@N array
    Of over/under-sampled frequencies 
    """
    lobe = round(scale * N / 2)
    delta = v_s / scale / N
    return np.array([i * delta for i in range(-lobe, lobe)])

#def fr(v_s, N):
#    return np.linspace(-v_s / 2, v_s / 2, num=N)    

# Primarily comes in handy for analysis 5.6

#def freq_range(v_s, N, W):
#    """
#    Return a N-length array
#    Where frequencies range between plus or minus W * v_s / 2
#    """
    #lobe = round(W * v_s / 2)
#    lobe = round(N / 2)
#    interval = W * v_s / N
#    return np.array([i * interval for i in range(-lobe, lobe)])

# this default line poises the script for lab 1 week 1 data collection
create_time(6.25e6)

def load_saves(filename):
    """
    Return a dictionary containing the arrays saved to
    the .npz file at @filename
    """
    a = np.load(filename, allow_pickle=True)
    return dict(zip(("{}".format(k) for k in a), (a[k] for k in a)))

def normalize(arr, actual_max):
    """
    Return a copy of @arr
    where each value has been scaled such that its max value
    corresponds to the given @actual_max
    """
    norm_c = actual_max / max(arr)
    return arr * norm_c

def plot(x, y, xnote, ynote, dual=False, noteA=None, noteB=None, logv=False):
    """
    General but highly volatile plotting function
    (because I do not like writing functions in the shell).    
    Mostly ad-hoc and designed for rapid plotting 
    of similar categories of data.
    """
    plt.figure(figsize=(3,3))
    plt.subplots_adjust(left=.2, bottom=.15, right=.95, top=.9)
    if dual:
        plt.plot(x, y[0], label=noteA)
        plt.plot(x, y[1], label=noteB)
        plt.legend()
    else:    
        plt.plot(x, y)   
    
    plt.xlabel(xnote, fontsize=12)
    plt.ylabel(ynote, fontsize=12)
    if logv:    
        plt.yscale('log')	

    plt.show()

def limited_plot(x, y, xnote, ynote, ycaps=None, logv=False):
    plt.figure(figsize=(3,3))
    plt.subplots_adjust(left=.2, bottom=.15, right=.95, top=.9)
    plt.plot(x, y)   
    
    plt.xlabel(xnote, fontsize=12)
    plt.ylabel(ynote, fontsize=12)
    if ycaps is not None:    
        plt.ylim(ycaps)
    if logv:    
        plt.yscale('log')	

    plt.show()

def plot_sample(arr, norm, pad=0, i=0, f=None):
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
    plt.plot(times[i + pad:f + pad], normalize(arr[i:f], norm))	    
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
