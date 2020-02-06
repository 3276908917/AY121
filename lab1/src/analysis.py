import ugradio
import numpy as np
import matplotlib.pyplot as plt

# I am assuming that we are continuing to use sample frequency 6.25e6 Hz
def power(sample, rate):
    """
    Plot power spectrum of @sample at @rate sampling rate
    and return Fourier transform as well as power array
    """
    f = ugradio.dft.dft(sample, vsamp=rate)
    P = np.abs(f[1]) ** 2
    plt.plot(f[0], P)
    plt.show()
    return f, P

def invf(transform):
    """
    Plot inverse Fourier transform of @transform
    and return the plotted array
    """
    i = ugradio.dft.idft(transform)
    plt.plot(i[0], i[1])
    plt.show()
    return i

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
