# todo: clean up invf and put the old version in review.py,
    # just like you did with power_plot and raw_power
# todo: expand power_plot to return values of multiple spikes
    # this will make it re-usable for the mixed signals section

import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

def hist_gauss(sample, b=64, correction=6):
    sigma = np.std(sample)
    count, bins, ignored = plt.hist(sample, b, density=True)
    plt.plot(bins, correction/(sigma * np.sqrt(2 * np.pi)) * np.exp(-bins ** 2 / (2 * sigma ** 2)), linewidth=2, color='r')
    plt.ylabel(r'Frequency')
    plt.xlabel('Voltage [?]') # normalization is not yet possible
    plt.show()

def freq_range(v_s, N, W=1):
    """
    Return a N-length array
    Where frequencies range between plus or minus W * v_s / 2
    """
    lobe = round(N / 2)
    interval = W * v_s / N
    return np.array([i * interval for i in range(-lobe, lobe)])

def power_barrage(bblock):
    power_book = []    
    #x = freq_range(srate, nsamps) / 10**6
    for c in bblock:
        f = np.fft.fft(c)
        P = np.abs(f) ** 2
        power_book.append(P)
    return power_book

# The normalization changes between trials?
# If so, that may explain the limited utility of naive averaging
def pp_skeleton(x, y, xBounds=None, yBounds=None, logv=False):
    fig = plt.figure()

    ax = fig.add_subplot(111)

    ax.plot(x, np.fft.fftshift(y))
    plt.xlabel('Frequency (MHz)')
    plt.ylabel(r'Magnitude-squared Voltage (V$^2$)')
    if xBounds is not None:
        plt.xlim(xBounds)
    if yBounds is not None:
        plt.ylim(yBounds)
    if logv:
        plt.yscale('log')
    plt.show()

def power_plot(sample, norm, srate=6.25e6, nsamps=16000, ifreq=None):
    """
    Clean, labeled plot of power spectrum of @sample
    corresponding to a signal at @ifreq (in MHz)
    with maximum amplitude @norm (in mV)
    sampled at rate @srate (default is our value for week 1 data collection)
    returns power spectrum and
        SINGLE frequency of largest power spike,
            so this function does not work as well for
            combined signals!
	Note: @ifreq is not used in any calculations,
		but rather is used to title the plot.
		Leave this as None, and the plot will not be titled.
    """
    fig = plt.figure()

    ax = fig.add_subplot(111)
    f = np.fft.fft(sample)
        #f = ugradio.dft.dft(sample, f=freqs, vsamp=srate)  
    P = np.abs(f) ** 2
    x = freq_range(srate, nsamps) / 10**6
    y = normalize(P, norm ** 2) / 1000 ** 2
    ax.plot(x, np.fft.fftshift(y))
    plt.xlabel('Frequency (MHz)')
    plt.ylabel(r'Magnitude-squared Voltage (V$^2$)')
    if ifreq is not None: # this is an awful way of giving a no-plot option
        plt.title('Power Spectrum: ' + str(ifreq) + ' MHz sinusoid')
    #plt.show()
    print(f[np.argmax(y)])
    return [x, y] #untested return... the x and y values
        # may not be aligned properly
