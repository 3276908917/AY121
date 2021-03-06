# from within a shell:
# execfile('plotter.py') for Python2
# exec(open('plotter.py').read()) for Python3


# git init
# git remote add origin master <link>


import ugradio
import numpy as np
import matplotlib.pyplot as plt

# abbreviation for "get"
def g(d):
	return ugradio.pico.capture_data('1V', divisor=d)

#dual get
def d():
	return ugradio.pico.capture_data('1V', divisor=2, dual_mode=True)

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

# np.savez('redos', re1=redo1, re2=redo2, re3=redo3)


