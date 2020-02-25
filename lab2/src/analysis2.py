# todo: clean up invf and put the old version in review.py,
    # just like you did with power_plot and raw_power
# todo: expand power_plot to return values of multiple spikes
    # this will make it re-usable for the mixed signals section

#import ugradio
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

def hist_gauss(sample, b=64, correction=6):
    sigma = np.std(sample)
    count, bins, ignored = plt.hist(sample, b, density=True)
    plt.plot(bins, correction/(sigma * np.sqrt(2 * np.pi)) * np.exp(-bins ** 2 / (2 * sigma ** 2)), linewidth=2, color='r')
    plt.ylabel(r'Frequency')
    plt.xlabel('Voltage [?]')
    plt.show()

def freq_range(v_s, N, W=1):
    """
    Return a N-length array
    Where frequencies range between plus or minus W * v_s / 2
    """
    #lobe = round(W * v_s / 2)
    lobe = round(N / 2)
    interval = W * v_s / N
    return np.array([i * interval for i in range(-lobe, lobe)])

# Hard-coded massive pack, based on Max's data
def power_barrage(rblock, iblock, srate=62.5e6, nsamps=15800):
    power_book = []    
    #x = freq_range(srate, nsamps) / 10**6

    for c in range(len(rblock)):
        for k in range(len(rblock[0])):
            sample = complex_combine(rblock[c][k], iblock[c][k])
            f = np.fft.fft(sample)
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

def invf(transform):
    """
    Plot inverse Fourier transform of @transform
    and return the plotted array

    NOT YET IMPLEMENTED with numpy
    """
    #i = ugradio.dft.idft(transform[1], transform[0])
    plt.plot(i[0], i[1])
    plt.show()
    return i

def zoom_comparison(sample, cut, freq, norm):
    """
    First plot: raw data (should be blocky and only of significance as a whole)
    Second plot: @cut-sample subsection of data, to make sure it fits expected shape
    Third plot: power spectrum. Eyeballing is not very helpful, so
        the function also returns the index associated with the
        power-maximizing frequency.

    NOT YET IMPLEMENTED with numpy
    """
    fig = plt.figure()
    ax1 = fig.add_subplot(311)

    ax1.plot(times[:len(sample)], normalize(sample, norm))
    plt.xlabel('Time (s)')
    plt.ylabel(r'Voltage (mV)')
    plt.title('Voltage vs Time: ' + str(freq) + ' MHz sinusoid')

    ax2 = fig.add_subplot(312)
    ax2.plot(times[:cut], normalize(sample[:cut], norm))
    plt.xlabel('Time (s)')
    plt.ylabel(r'Voltage (mV)')
    plt.title('First ' + str(cut) + ' samples: ' + str(freq) + ' MHz sinusoid')

    ax3 = fig.add_subplot(313)
    return 1 / 0
    f = ugradio.dft.dft(sample, vsamp=6.25e6)
    P = np.abs(f[1]) ** 2
    ax3.plot(f[0], normalize(P, norm ** 2))
    plt.xlabel('Frequency (Hz)')
    plt.ylabel(r'Power (nV$^2$)')
    plt.title('Power Spectrum: ' + str(freq) + ' MHz sinusoid')

    plt.show()
    return f[0][np.argmax(P)]    

# Voltage spectrum
def volt_spec(sample, ifreq, norm, srate=6.25e6):
    """
    Plot 1: real and imaginary components of voltage spectrum of sample.
    Plot 2: just the real part
    Plot 3: just the imaginary part
    @ifreq input frequency
    @norm maximum observed amplitude of input signal
    @srate sampling frequency, default is lab 1 week 1
    
    Unless there is great overlap between imaginary and real components,
        use mini_volt_spec instead

    NOT YET IMPLEMENTED with numpy
    """
    
    fig = plt.figure()
    ax1 = fig.add_subplot(311)

    return 1 / 0
    f = ugradio.dft.dft(sample, vsamp=srate)

    ax1.plot(f[0], normalize(np.real(f[1]), norm), label='real component')
    ax1.plot(f[0], normalize(np.imag(f[1]), norm), label='imaginary component')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel(r'Voltage (mV)')
    plt.title('Voltage Spectrum: ' + str(ifreq) + ' MHz sinusoid')

    ax1.legend(bbox_to_anchor=(1, 1))

    ax2 = fig.add_subplot(312)

    ax2.plot(f[0], normalize(np.real(f[1]), norm))
    plt.xlabel('Frequency (Hz)')
    plt.ylabel(r'Voltage (mV)')
    plt.title('Real Voltage Spectrum: ' + str(ifreq) + ' MHz sinusoid')

    ax3 = fig.add_subplot(313)

    ax3.plot(f[0], normalize(np.imag(f[1]), norm))
    plt.xlabel('Frequency (Hz)')
    plt.ylabel(r'Voltage (mV)')
    plt.title('Imaginary Voltage Spectrum: ' + str(ifreq) + ' MHz sinusoid')

    plt.show()

def mini_volt_spec(sample, norm, srate=6.25e6, ifreq=None):
    """
    Simultaneously plot real and imaginary components of voltage spectrum of sample.
    @ifreq input frequency
    @norm maximum observed amplitude of input signal
    @srate sampling frequency, default is lab 1 week 1
    Returns the Fourier transform.
    If there is excessive overlap between imaginary and real components,
        use volt_spec instead

    NOT YET IMPLEMENTED with numpy
    """
    
    fig = plt.figure()
    ax1 = fig.add_subplot(111)

    return 1 / 0
    f = ugradio.dft.dft(sample, vsamp=srate)    

    x = f[0] / 10 ** 6
    re = normalize(np.real(f[1]), norm)
    im = normalize(np.imag(f[1]), norm)

    ax1.plot(x, re, label='real component')
    ax1.plot(x, im, label='imaginary component')
    plt.xlabel('Frequency (MHz)')
    plt.ylabel(r'Voltage (mV)')
    if ifreq is not None: # this is an awful way of giving a no-plot option
        plt.title('Voltage Spectrum: ' + str(ifreq) + ' MHz sinusoid')

    ax1.legend(bbox_to_anchor=(1, 1))

    plt.show()
    return [x, [re, im]]

def full_comparison(sample, rate):
    """
    Plot 1: original array @sample vs time 'times'
    Plot 2: power spectrum of @sample (frequency vs square volt)
    Plot 3: inverse fourier transform of power spectrum

    NOT YET IMPLEMENTED with numpy
    """
    
    fig = plt.figure()
    ax1 = fig.add_subplot(311)
    ax1.plot(times[:len(sample)], sample)

    ax2 = fig.add_subplot(312)
    return 1 / 0
    f = ugradio.dft.dft(sample, vsamp=rate)
    P = np.abs(f[1]) ** 2
    ax2.plot(f[0], P)

    # It is not clear to me that this next section is actually beneficial;
    # When taking the power spectrum we are bound to lose information
    # in the real and imaginary coefficients.

    ax3 = fig.add_subplot(313)
    i = ugradio.dft.idft(P)

    shift = np.append(i[1][8000:16000], i[1][0:8000])

    print(len(shift), len(sample))

    ax3.plot(times[:len(sample)], shift)
    plt.show()

    # In conclusion, I think that the manual is suggesting
    # that we take the inverse Fourier transform directly
    # of the Fourier transform (no extra math in between).