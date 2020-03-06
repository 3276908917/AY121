# todo: clean up invf and put the old version in review.py,
    # just like you did with power_plot and raw_power
# todo: expand power_plot to return values of multiple spikes
    # this will make it re-usable for the mixed signals section

import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

def hist_gauss(sample, b=64, correction=6):
    '''
    Given an array of values @sample,
    return a histogram consisting of @b bins (64 by default)
    and over-plot what a zero-mean Gaussian would look like with the same standard deviation,
    and with an amplitude correction factor given by @correction (6 by default). 
    '''
    sigma = np.std(sample)
    count, bins, ignored = plt.hist(sample, b, density=True)
    plt.plot(bins, correction/(sigma * np.sqrt(2 * np.pi)) * np.exp(-bins ** 2 / (2 * sigma ** 2)), linewidth=2, color='r')
    plt.ylabel(r'Frequency')
    plt.xlabel('Voltage [?]') # normalization is not yet possible
    plt.show()

def power_barrage(bblock):
    '''
    Given a big block of arrays @bblock (each array is a data capture),
    return a big block of their corresponding power spectra.
    '''
    power_book = []    
    for c in bblock:
        f = np.fft.fft(c)
        P = np.abs(f) ** 2
        power_book.append(P)
    return power_book

# The normalization changes between trials?
# If so, that may explain the limited utility of naive averaging
def pp_skeleton(x, y, xBounds=None, yBounds=None, logv=False, xLabel='Frequency (MHz)',
    yLabel = r'Magnitude-squared Voltage (V$^2$)'):
    '''
    Plot @y versus @x where @y is shifted to center the zero frequency
        (more at numpy.fft.fftshift)
    The plot is automatically labeled such that the
        x-axis corresponds to frequency in megahertz
        y-axis corrosponds to magnitude-squared voltage in square volts
    One can include tuples @xBounds and @yBounds to zoom the graph appropriately
        where each tuple is in the form (lower_bound, upper_bound).
    One can ahead-of-time specify the y-axis to be logarithmic
        with @logv=True. However, it is not clear whether this functionality is useful,
        because one can simply press L (while focused on a matplotlib plot)
        to achieve the same effect.
    '''
    fig = plt.figure(figsize=(6,3))
    plt.subplots_adjust(left=.15, bottom=.15, right=.95, top=.9)

    ax = fig.add_subplot(111)

    ax.tick_params(axis="x", labelsize=12)
    ax.tick_params(axis="y", labelsize=12)

    ax.plot(x, np.fft.fftshift(y))
    plt.xlabel(xLabel, fontsize=12)
    plt.ylabel(yLabel, fontsize=12)
    if xBounds is not None:
        plt.xlim(xBounds)
    if yBounds is not None:
        plt.ylim(yBounds)
    if logv:
        plt.yscale('log')
    plt.show()

# Double power plot. We share an axis and stack vertically
def double_pp(x, y1, y2, xBounds=None, yBounds=None, logv=False):
    fig, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=True)

    #fig.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axes
    #plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    #plt.grid(False)

    plt.xlabel('Frequency (MHz)', fontsize=12)
    fig.text(0, 0.5, r'Magnitude-squared Voltage (V$^2$)',
             va='center', rotation='vertical', fontsize=12)

    ax1.tick_params(axis="x", labelsize=12)
    ax2.tick_params(axis="x", labelsize=12)
    ax1.tick_params(axis="y", labelsize=12)
    ax2.tick_params(axis="y", labelsize=12)
    

    ax1.plot(x, np.fft.fftshift(y1))
    ax2.plot(x, np.fft.fftshift(y2))    
    
    if xBounds is not None:
        plt.xlim(xBounds)
    if yBounds is not None:
        plt.ylim(yBounds)
    if logv:
        plt.yscale('log')
        
# over-plotting of two distributions on the same graph.
def over_pp(x, y1, y2, y1L, y2L, xBounds=None, yBounds=None, logv=False):
    fig = plt.figure(figsize=(6,3))
    plt.subplots_adjust(left=.15, bottom=.15, right=.95, top=.9)

    ax = fig.add_subplot(111)

    ax.tick_params(axis="x", labelsize=12)
    ax.tick_params(axis="y", labelsize=12)

    ax.plot(x, np.fft.fftshift(y1), label=y1L)
    ax.plot(x, np.fft.fftshift(y2), label=y2L)
    plt.xlabel('Frequency (MHz)', fontsize=12)
    plt.ylabel(r'Magnitude-squared Voltage (V$^2$)', fontsize=12)
    if xBounds is not None:
        plt.xlim(xBounds)
    if yBounds is not None:
        plt.ylim(yBounds)
    if logv:
        plt.yscale('log')

    ax.legend(bbox_to_anchor=(1, 1))
        
    plt.show()

# each gauss should be a triple (amp, avg, sig)
def tau_cannon(x, y, gauss1, gauss2, xBounds=None, yBounds=None):
        fig = plt.figure(figsize=(6,3))
    plt.subplots_adjust(left=.15, bottom=.15, right=.95, top=.9)

    ax = fig.add_subplot(111)

    ax.tick_params(axis="x", labelsize=12)
    ax.tick_params(axis="y", labelsize=12)

    ax.plot(x, np.fft.fftshift(y1), label=y1L)
    ax.plot(x, np.fft.fftshift(y2), label=y2L)
    plt.xlabel('Frequency (MHz)', fontsize=12)
    plt.ylabel(r'Magnitude-squared Voltage (V$^2$)', fontsize=12)
    if xBounds is not None:
        plt.xlim(xBounds)
    if yBounds is not None:
        plt.ylim(yBounds)
    if logv:
        plt.yscale('log')

    ax.legend(bbox_to_anchor=(1, 1))
        
    plt.show()

def power_plot(sample, norm, srate=6.25e6, nsamps=16000, ifreq=None):
    '''
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
    '''
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    f = np.fft.fft(sample)
    P = np.abs(f) ** 2
    
    x = freq_range(srate, nsamps) / 10**6
    y = normalize(P, norm ** 2) / 1000 ** 2
    
    ax.plot(x, np.fft.fftshift(y))
    plt.xlabel('Frequency (MHz)')
    plt.ylabel(r'Magnitude-squared Voltage (V$^2$)')
    
    if ifreq is not None: # this is an awful way of giving a no-plot option
        plt.title('Power Spectrum: ' + str(ifreq) + ' MHz sinusoid')

    plt.show()
    print(f[np.argmax(y)])
    return [x, y] # untested return...
        # ...the x and y values may not be aligned properly
