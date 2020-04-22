import ugradio
import numpy as np
import matplotlib.pyplot as plt
from ugradio import timing

def fourier_skeleton(x, y, xBounds=None, yBounds=None, logv=False,
    xLabel='Frequency (Hz)', yLabel = r'Voltage (mV)', rude_filter = False):
    fig = plt.figure(figsize=(6,3))
    plt.subplots_adjust(left=.15, bottom=.15, right=.95, top=.9)

    ax = fig.add_subplot(111)

    ax.tick_params(axis="x", labelsize=12)
    ax.tick_params(axis="y", labelsize=12)

    if rude_filter:
        y[0:1] = 0

    ax.plot(x, np.real(np.fft.fftshift(y * 1000)), label='real')
    ax.plot(x, np.imag(np.fft.fftshift(y * 1000)), label='imaginary')
    
    plt.xlabel(xLabel, fontsize=12)
    plt.ylabel(yLabel, fontsize=12)
    if xBounds is not None:
        plt.xlim(xBounds)
    if yBounds is not None:
        plt.ylim(yBounds)
    if logv:
        plt.yscale('log')
    plt.show()

    # Now we want to find the locations of the peaks
    print(np.abs(x[np.argmax(np.fft.fftshift(y))]))

# hard-code city
def collect_fringes(volts):
    fringes = []
    # We evaluate by ten minute segments (10 * 60 = 600 seconds per)
    x = fr(1, 600)
    P = lambda r : np.abs(np.fft.fft(r)) ** 2
    for i in range(600, len(volts), 600):
        y = volts[i - 600:i]
        y = P(y)
        y[0:1] = 0 # naive filtering approach
        y = np.fft.fftshift(y)
        fringes.append(np.abs(x[np.argmax(y)]))
    #print(fringes)
    return fringes

def pp3(x, y, xBounds=None, yBounds=None, logv=False, xLabel='Frequency (MHz)',
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

'''
Notes to self:
    sunburst seems to stop at about 1.5853214e9
        and resumes at about        1.58535819e9
    indices 11k to 42k seem safe

    since we took data at 1 Hz, the recommended 10 minutes
        will correspond to 600 indices.

'''
