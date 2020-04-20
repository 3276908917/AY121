import ugradio
import scipy.signal 
from astropy.io import fits
import numpy as np
import ugradio.coord
from scipy import fft
import scipy.stats as stats
import scipy
import matplotlib.pyplot as plt

def mean(signal):
	return np.mean(signal)

def root_mean_square(signal):
	return np.sqrt(mean(signal)**2)

def standard_deviation(signal):
	return np.std(signal)

def power_spectrum(signal):
	fourier_transform= np.fft.fft(signal)
	return np.abs(fourier_transform)**2

def variance(signal):
	return np.var(signal)

def fourier(signal):
	return np.fft.fft(signal)
	
def load_saves(filename):
	a = np.load(filename)
	return dict(zip(("{}".format(k) for k in a), (a[k] for k in a)))





def Read_leushner_header(file_name):
	
	""" Get data from file name """
	
	leushner_data = fits.open(str(file_name))
	
	return leushner_data


def Data_header(file_name):
	
	""" Get header from file name """
	
	leushner_header = fits.open(str(file_name))
	
	return leushner_header[0].header
	

	
def Plot_leushner_raw_data(file_name):
	
	"""
	Graph animmated plots of both auto0_real and auto0_real over all the samples N
	"""
	f = fits.open(str(file_name))
	
	N = f[0].header['NSPEC']
	
	for p in np.arange(1, N):
		
		plt.ion()
		
	
		plt.plot(f[p].data['auto0_real'], label=' auto0_real')
		plt.plot(f[p].data['auto1_real'], label= 'auto1_real')
	
		plt.xlabel('samples')
		plt.ylabel('Amplitude ($V^2$)')
	#	plt.title(' Power spectrum ')
		plt.legend()
		
		plt.draw()
		plt.pause(1)


		
			
	plt.ioff()
					
	return f


def Average_spectrum(leushner_data, N=10, polarization_number='first'):
	
	"""
	Get Average of signals based on numbers of spectras N and polarization number first/second 
		
	"""
	
	average_spectrum = []


	if polarization_number == "first":
		
		name = 'auto0_real'
		
		for i in np.arange(1, N):
#			print ('leushner_data[i].data[name]', leushner_data[i].data[name])
			spectrum_data = leushner_data[i].data[name]
			average_spectrum.append( spectrum_data )
			
			
		Average = np.mean(average_spectrum, axis = 0)
			
		return Average
			
				
			
	elif polarization_number == "second":
		
		name = 'auto1_real'
		
		for i in np.arange(1, N):
			
			spectrum_data = leushner_data[i].data[name]
			average_spectrum.append(spectrum_data )
			
			
		Average = np.mean(average_spectrum, axis = 0)
			
		return Average

def plot_average(file_name, polarization_number= 'first' ):
	
	
	leushner_data = Read_leushner_header(file_name)
			
	header = Data_header(file_name)

	N = header['NSPEC']
		
		
	#average of data, either first or second 
	average_spectra = Average_spectrum(leushner_data, N, polarization_number)
	
	
	print ('Average_spectra',average_spectra)
	
	#frequency responce of our AVERAGE data spectrum
	freq = Frequency_spectrum(header, average_spectra)
	
	plt.plot(freq, average_spectra, label=file_name)
	plt.xlabel('Frequency')
	plt.ylabel('$V^2$')
	plt.legend()
#	plt.show()

				
def Frequency_spectrum(header, spectrum_data):
	
	"""
	Changes samples axis of the spectrum data into frequency axis	
	"""
		
	vsamp = header["SAMPRATE"]
		
	t = 1/vsamp
	
	freq = (np.fft.fftshift(np.fft.fftfreq(len(spectrum_data), t ))) / 1e6 + 1420  # where 1420 is the total of the LO frequencies in MHz that shifted  down our HI signal, I added this to re-shift the received signal toward frequency 1420.4 mhz

	return freq


def base_fit_average(average_spectra, frequency):
	
	"""	
	Filter our average data, then use 1d interpolation on the data within boundaries 
	
	"""
		
	filtered_spectra = scipy.signal.medfilt(average_spectra, kernel_size=9) 
	max_index = np.argmax(filtered_spectra)
	
	bandwidth = 200
	mask = np.ones(len(filtered_spectra), dtype=bool)
	
	for i in range(len(mask)):
		if  (max_index - bandwidth) < i < (max_index + bandwidth) :
			mask[i] = False
			
	baseline_interpolation = scipy.interpolate.interp1d(frequency[mask], filtered_spectra[mask])
	
	
#	plt.plot(frequency,baseline_interpolation(frequency))	
#	plt.plot(frequency, filtered_spectra, c='b')
#	plt.plot(frequency, average_spectra, marker='x', c='r')
#	plt.show()
	

	return baseline_interpolation(frequency)


def doppler_velocity_correction(header, frequency, velocity_correction='off'):
	
	"""
	calculates the correction needed to shift the doppler velocity to the LSR frame. LSR frame Local standard of rest frame, which follows the mean motion of material in the Milky Way in the neighborhood of the Sun.
	
	We use the header from the collected data each time to get the ra, dec and jd; and use them to get the doppler correction, then we get the doppler velocity	
	"""
	#parameters
	HI_freq = 1420.405751  # frequency of our HI signal in Mhz
	c = 3e8  #speed of light in km/s 
		
		
	#Get doppler correction
	ra, dec, jd = header["RA"], header["DEC"], header["JD"]
	print ('ra, dec, jd = ',ra, dec, jd)
		
	

#	calculating doppler velocity
	delta_freq = (frequency) - HI_freq
	doppler_velocity = delta_freq / HI_freq * c
	
	

	#turn on Doppler correction
	if velocity_correction=='on':
				
		Doppler_correction = ugradio.doppler.get_projected_velocity(ra, dec, jd) 

		print('Doppler_correction',Doppler_correction)
		Velocity =  np.array(doppler_velocity) - Doppler_correction.value
	
	#turn off Doppler correction
	elif velocity_correction=='off':
		Velocity =  np.array(doppler_velocity)
		 
	return  Velocity



def Gain(signal_noise_on, signal_noise_off, number_spectra, Noise_temp, Sky_temp):
	
	"""
	Calculate the Gain (G) like in Lab 2, uses the signal with and without noise, temperature of the sky and temprature from the noise in kelvin
	"""
	
	
	noise_on_spectra = Average_spectrum(signal_noise_on, number_spectra )
	
	noise_off_spectra = Average_spectrum(signal_noise_off, number_spectra )
	
	

#	Noise_temp - Sky_temp = 300 #kelvin

	gain = ((Noise_temp - Sky_temp) / np.sum(noise_on_spectra - noise_off_spectra)) * np.sum(noise_off_spectra)
	
	print('Gain',gain)
	
	return gain



	

def final_calibration(file_name, Polarization_number='first', Gain=1, velocity_correction='off' ):
	
	
	"""
	-Calculate the gain first Gain(noise_on, noise_off, Noise_temp, Sky_temp) then use the final value of the gain for the calibration
	
	-Calibrate the leushner data, change it from frequency Mhz vs amplitude V^2 to temperature spectra (kelvin) with respect to doppler velocity m/s (LSR frame)
	
	"""
	
	gain = 300 # will calculate it later once
	
	leushner_data = Read_leushner_header(file_name)
		
	header = Data_header(file_name)

	N = header['NSPEC']
	
	
	#average of data, either first or second 
	average_spectra = Average_spectrum(leushner_data, N, Polarization_number)
	print ('Average_spectra', average_spectra)
	#frequency responce of our AVERAGE data spectrum
	freq = Frequency_spectrum(header, average_spectra)

	Temperature = (average_spectra - base_fit_average(average_spectra, freq)) * Gain
	print ('Temperature', Temperature)
	
	doppler_velocity = doppler_velocity_correction(header, freq, velocity_correction)
	
	
	return Temperature, doppler_velocity
	
	

fon= fits.open('thirty_plane_on.fits')
N = fon[0].header['NSPEC']
	

foff= fits.open('thirty_plane_off.fits')
N = foff[0].header['NSPEC']


gain = Gain(fon, foff, N, 301, 1)
#calculate the gain:
print ('Gain', Gain(fon, foff, N, 301, 1))


#comparing on/off data, by inputing file name and polarization number 

plot_average('thirty_plane_on.fits')
plot_average('thirty_plane_off.fits')
plt.show()

plot_average('thirty_plane_on.fits','second')
plot_average('thirty_plane_off.fits','second')
plt.show()


#plot all rawa data from file
#Plot_leushner_raw_data('thirty_plane_on.fits')





#calibration, (ugradio.doppler.get_projected_velocity(ra, dec, jd)) outputs the doppler velocity and the unit, this format blocks the velocity correction calculation when it's on.

data = final_calibration('thirty_plane_on.fits', velocity_correction='off', Gain= -1*gain)
plt.plot(data[1], data[0])
plt.show()
