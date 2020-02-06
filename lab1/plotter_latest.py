# from within a shell:
# execfile('plotter.py')

import ugradio
import numpy as np
import matplotlib.pyplot as plt

# abbreviation for "get"
def g():
        """
        Abbreviation function.
        Acquire data through pico sampler. 1V range, divisor=10.
        """
	return ugradio.pico.capture_data('1V', divisor=10)

times = []

for i in range(0, 16000):
	times.append(i / 6.25e6) # time intervals given by sampling frequency

def plot(arr):
	plt.plot(times[:len(arr)], arr)
	plt.show()

def collect_noise():
	noises = []
	for c in range (0, 32):
		noises.append(g())
	return noises

# returns dictionary with arrays of collected data
def load_saves(filename):
	a = np.load(filename)
	return dict(zip(("{}".format(k) for k in a), (a[k] for k in a)))

# I am assuming that we are continuing to use sample frequency 6.25e6 Hz
def power(sample):
	f = ugradio.dft.dft(sample, vsamp=6.25e6) # Is it discouraged to avoid the optional args?	
	P = np.abs(f[1]) ** 2
	plt.plot(f[0], P)
	plt.show()
	return P

def power(sample):
	f = ugradio.dft.dft(sample, vsamp=6.25e6)
	P = np.abs(f[1]) ** 2
	return f, P

def invf(spectrum):
	i = ugradio.dft.idft(spectrum)
	plt.plot(i[0], i[1])
	plt.show()
	return i

# I need a data normalization function:
# normalize(arr, max)
        # calculates the ratio between the max value of arr and the given param @max
        # It then returns a copy of arr normalized by this ratio

def normalize(arr, actual_max):
        norm_c = actual_max / max(arr)
        return arr * norm_c
        
"""
First plot: raw data (should be blocky and only of significance as a whole)
Second plot: @cut-sample subsection of data, to make sure it is sinusoidal
Third plot: power spectrum. Eyeballing is not very helpful, so
        the function also returns the index associated with the
        power-maximizing frequency.
"""
def zoom_comparison(sample, cut, freq, norm):
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

"""
Simultaneously plot real and imaginary components of voltage spectrum of sample.
"""
def voltage_spectrum(sample, freq, norm):
	fig = plt.figure()
	ax1 = fig.add_subplot(311)

	f = ugradio.dft.dft(sample, vsamp=6.25e6)

	ax1.plot(f[0], normalize(np.real(f[1]), norm), label='real component')
	ax1.plot(f[0], normalize(np.imag(f[1]), norm), label='imaginary component')
	plt.xlabel('Frequency (Hz)')
	plt.ylabel(r'Voltage (mV)')
	plt.title('Voltage Spectrum: ' + str(freq) + ' MHz sinusoid')

	ax1.legend(bbox_to_anchor=(1, 1))

	ax2 = fig.add_subplot(312)

	ax2.plot(f[0], normalize(np.real(f[1]), norm))
	plt.xlabel('Frequency (Hz)')
	plt.ylabel(r'Voltage (mV)')
	plt.title('Real Voltage Spectrum: ' + str(freq) + ' MHz sinusoid')

	ax3 = fig.add_subplot(313)

	ax3.plot(f[0], normalize(np.imag(f[1]), norm))
	plt.xlabel('Frequency (Hz)')
	plt.ylabel(r'Voltage (mV)')
	plt.title('Imaginary Voltage Spectrum: ' + str(freq) + ' MHz sinusoid')

	plt.show()

def mini_voltage_spectrum(sample, freq, norm):
	fig = plt.figure()
	ax1 = fig.add_subplot(111)

	f = ugradio.dft.dft(sample, vsamp=6.25e6)

	ax1.plot(f[0], normalize(np.real(f[1]), norm), label='real component')
	ax1.plot(f[0], normalize(np.imag(f[1]), norm), label='imaginary component')
	plt.xlabel('Frequency (Hz)')
	plt.ylabel(r'Voltage (mV)')
	plt.title('Voltage Spectrum: ' + str(freq) + ' MHz sinusoid')

	ax1.legend(bbox_to_anchor=(1, 1))

	plt.show()

def mini_power_spectrum(sample, freq, norm):
	fig = plt.figure()

	ax3 = fig.add_subplot(111)
	f = ugradio.dft.dft(sample, vsamp=6.25e6)
	P = np.abs(f[1]) ** 2
	ax3.plot(f[0], normalize(P, norm ** 2))
	plt.xlabel('Frequency (Hz)')
	plt.ylabel(r'Power (nV$^2$)')
	plt.title('Power Spectrum: ' + str(freq) + ' MHz sinusoid')

	plt.show()
	
	return f[0][np.argmax(P)]

def full_comparison(sample):
	fig = plt.figure()
	ax1 = fig.add_subplot(311)
	ax1.plot(times[:len(sample)], sample)

	ax2 = fig.add_subplot(312)
	f = ugradio.dft.dft(sample, vsamp=6.25e6)
	P = np.abs(f[1]) ** 2
	ax2.plot(f[0], P)

	ax3 = fig.add_subplot(313)
	i = ugradio.dft.idft(P)

	shift = np.append(i[1][8000:16000], i[1][0:8000])

	print(len(shift), len(sample))

	ax3.plot(times[:len(sample)], shift)
	plt.show()

# Here is an example of a labeled save: 
# np.savez('bundle1', trial1=t1, trial2=t2, trial3=t3, trial4=t4, trial5=t5, trial6_1=t6_1, trial6_2=t6_2, trial6_3=t6_3, trial7=t7, noises=n) 


