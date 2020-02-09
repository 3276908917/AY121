# todo: clean up invf and put the old version in review.py,
    # just like you did with power_plot and raw_power
# todo: expand power_plot to return values of multiple spikes
    # this will make it re-usable for the mixed signals section

import ugradio
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

def power_plot(sample, norm, srate=6.25e6, ifreq=None):
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
    f = ugradio.dft.dft(sample, vsamp=srate)
    P = np.abs(f[1]) ** 2
    x = f[0] / 10 ** 6
    y = normalize(P, norm ** 2) / 1000 ** 2
    ax.plot(x, y)
    plt.xlabel('Frequency (MHz)')
    plt.ylabel(r'Power (V$^2$)')
    if ifreq is not None: # this is an awful way of giving a no-plot option
        plt.title('Power Spectrum: ' + str(ifreq) + ' MHz sinusoid')

    plt.show()
    print(f[0][np.argmax(y)])
    return [x, y]

def invf(transform):
    """
    Plot inverse Fourier transform of @transform
    and return the plotted array
    """
    i = ugradio.dft.idft(transform) # remember to use f[1], not f, and not f[0]
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
    """
    
    fig = plt.figure()
    ax1 = fig.add_subplot(311)

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
    """
    
    fig = plt.figure()
    ax1 = fig.add_subplot(111)

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
    """
    
    fig = plt.figure()
    ax1 = fig.add_subplot(311)
    ax1.plot(times[:len(sample)], sample)

    ax2 = fig.add_subplot(312)
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
