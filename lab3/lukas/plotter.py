import numpy as np
import matplotlib.pyplot as plt

## Plotting functions removed from the object

def fourier_frequency(signal, sampling_frequency):
    '''
    TODO: this function does NOT plot the real
    and imaginary components of the spectrum,
    so the current output is a cast!
    '''
    signal_fft = np.fft.fft(signal)
    time_step = 1 / sampling_frequency
    freqs = np.fft.fftfreq(signal_fft.size, time_step)
    index = np.argsort(freqs)

    plt.plot(freqs[index] / 1e6, ps[index])
    plt.xlabel('Frequency (MHz)')
    plt.ylabel('Amplitude (V)')
    plt.title('Fourier transform')

    plt.legend()
    plt.show()

    return freqs[index], signal_fft[index]

def power_spectrum(signal):
    fourier_transform = np.fft.fft(signal)
    return np.abs(fourier_transform) ** 2

def spectrum_frequency(signal, sampling_frequency):
    ps = power_spectrum(signal)
    time_step = 1 / sampling_frequency
    freqs = np.fft.fftfreq(signal.size, time_step)
    index = np.argsort(freqs)

    #plot the frequency power spectrum
    plt.plot(freqs[index] / 1e6, ps[index])
    plt.xlabel('Frequency (MHz)')
    plt.ylabel(r'Amplitude ($V^2$)')
    plt.title('Power spectrum')
    # plt.yscale('log')

    plt.legend()
    plt.show()

    return freqs[index], ps[index]

## end section
