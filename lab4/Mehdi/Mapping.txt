import ugradio
import scipy.signal 
from astropy.io import fits
import numpy as np
import ugradio.coord
from scipy import fft
import scipy.stats as stats
import scipy
import matplotlib.pyplot as plt
import bisect
from mpl_toolkits.basemap import Basemap
from scipy.signal import find_peaks



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

def convert_kpc_to_km(kpc):
	return kpc * 3.086e16

def convert_km_to_kpc(km):
	return km / 3.086e16

def Magnitude(x, y):
	return np.sqrt(x**2 + y**2)




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
	HI_freq = 1420.40575  # frequency of our HI signal in Mhz
	c = 3e5  #speed of light in km/s 
		
		
	#Get right ascension, declination and julian time
	ra, dec, jd = header["RA"], header["DEC"], header["JD"]
	print ('ra, dec, jd = ',ra, dec, jd)
		
	

#	calculating doppler velocity
	delta_freq = (frequency) - HI_freq
	doppler_velocity = delta_freq / HI_freq * c
	
	

	#turn on Doppler correction
	if velocity_correction=='on':
				
		Doppler_correction = ugradio.doppler.get_projected_velocity(ra, dec, jd) 

		print('Doppler_correction',Doppler_correction)
		print('doppler_velocity',doppler_velocity)
		
		Velocity =  np.array(doppler_velocity) - (Doppler_correction.value/1e3) #convertion Doppler correction from m/s to km/s
	
	#turn off Doppler correction
	elif velocity_correction=='off':
		Velocity =  np.array(doppler_velocity)
		 
	return  Velocity



def Gain(signal_noise_on, signal_noise_off, number_spectra=10, Noise_temp=301, Sky_temp=1):
	
	"""
	Calculate the Gain (G) like in Lab 2, uses the signal with and without noise, temperature of the sky and temprature from the noise in kelvin
	"""
	
	
	noise_on_spectra = Average_spectrum(signal_noise_on, number_spectra )
	
	noise_off_spectra = Average_spectrum(signal_noise_off, number_spectra )
	
#	Noise_temp - Sky_temp = 300 #kelvin

	gain = ((Noise_temp - Sky_temp) / np.sum(noise_on_spectra - noise_off_spectra)) * np.sum(noise_off_spectra)
		
	return gain



	

def final_calibration(file_name, Polarization_number='first', Gain=1, velocity_correction='off' ):
		
	"""
	- Calculate the gain first Gain(noise_on, noise_off, Noise_temp, Sky_temp) then use the final value of the gain for the calibration
	
	- Calibrate the leushner data, change it from frequency Mhz vs amplitude V^2 to temperature spectra (kelvin) with respect to doppler velocity m/s (LSR frame)
	
	"""
	
#	gain = 300 # will calculate it later once
	
	
	
	
	leushner_data = Read_leushner_header(file_name)
		
	header = Data_header(file_name)

	N = header['NSPEC']
	
	
	#average of data, either first or second 
	average_spectra = Average_spectrum(leushner_data, N, Polarization_number)
	
	print ('Average_spectra', average_spectra)
	#frequency responce of our AVERAGE data spectrum
	freq = Frequency_spectrum(header, average_spectra)

	Temperature = 100* (average_spectra - base_fit_average(average_spectra, freq) ) * Gain
	print ('Temperature', Temperature)
	
	doppler_velocity = doppler_velocity_correction(header, freq, velocity_correction)
	
	
	return doppler_velocity, Temperature 
		
	


def plot_all_attempts(lo_freq=635, file_name='attempt_all', x_size=360 , y_size=180 ,velocity_correction='on', Gain= 800 ):
		
	degree_spacing = 2/2

# Galactic longtitude and latitude
	lontitude_min= -9
	longtitude_max = 250
	Longtitude_Spacing = (longtitude_max - (lontitude_min)) / degree_spacing
	
	longitude_range = np.linspace(lontitude_min, longtitude_max, int(Longtitude_Spacing)+1)
	 
	print ('longitude_range',len(longitude_range))
	Grid = np.zeros(shape=(180,360))
	Grid_velocity = np.zeros(shape=(180,360))
		
	for l in longitude_range:
		plt.ion()		
		try:
			
			print('Data cycle auto')						
			data_noisy_name = 'cycle_auto_' + str(float(int(l))) +'_degrees_'+ str(lo_freq) +'MHz_noisy.fits'	
			data_quiet_name = 'cycle_auto_' + str(float(int(l))) +'_degrees_'+ str(lo_freq) +'MHz_quiet.fits'		
			
			print (l )			
			print (data_noisy_name)	
			print (data_quiet_name)
			
			data_on = Read_leushner_header(str(data_noisy_name))
			data_off = Read_leushner_header(str(data_quiet_name))
			H = Data_header(str(data_quiet_name))

			
			gain = 300 / np.sum(Average_spectrum(data_on, N=10, polarization_number='first') - Average_spectrum(data_off, N=10, polarization_number='first')) * np.sum(Average_spectrum(data_off, N=10, polarization_number='first'))
			
#			print (Average_spectrum(data_on, N=10, polarization_number='first'))
#			print (Average_spectrum(data_off, N=10, polarization_number='first'))

			
			ra, dec = H["RA"], H["DEC"]					
			print ('ra = ', ra, ' dec = ', dec)	

			N = H['NSPEC'] #number spectra
			print (N)
		
			
			doppler_velocity, Temperature  = final_calibration(data_quiet_name, Polarization_number='first', Gain= gain, velocity_correction='on' )
			
			
			plt.plot(doppler_velocity[2000:5000], Temperature[2000:5000] )
#			plt.plot(doppler_velocity[2000:5000])
#			plt.plot(Temperature[2000:5000] )
		
			
			
						
#			peaks, _ = find_peaks(doppler_velocity, height=0)
#		
#			plt.plot(peaks, doppler_velocity[peaks])
#			plt.show()
			
#			print(max(Temperature[peaks]))

			index = np.where(Temperature == max(Temperature[0:5500]))
			
#			print(index_max)
			print ('HI_temperature = ', max(Temperature), np.where(Temperature == max(Temperature)),  'HI_Velocity = ', doppler_velocity[index[0]][0])
			#Plot average plots of all data in the file within the longtitude we collected
#			
			#Plot Temperatude vs Doppler velocity of all data in the file within the longtitude we collected

			
			m = Basemap(projection='cyl') # Cylindrical projection
		
			Ra, Dec = m.makegrid(x_size , y_size)
			
			Ra_list = np.arange(0,360,10) # draw longitude lines every 30 degrees
			#lat_lines = np.arange(14,85,30) # draw latitude lines every 30 degrees
			dec_list = np.arange(-90,90,10)
			
#			print ('m', m.value())
			
			Grid[int(dec)][int(ra)]= max(Temperature[0:5500])
			Grid_velocity[int(dec)][int(ra)]= doppler_velocity[index[0]][0]
			
			

			
#			plt.imshow(Grid, cmap='plasma', interpolation='gaussian')
#			#	.drawmeridians(Ra_list, labels=[0,0,0,1], labelstyle='+/-')
#			#	m.drawparallels(dec_list, labels=[0,1,0,0], labelstyle='+/-')
#			plt.colorbar()
			plt.draw()
			plt.pause(.1)
			plt.close()

			
		
						
			
		except:
			
			print('Data cycle 2')
			data_noisy_name = 'cycle2_' + str(float(int(l))) +'_degrees_'+ str(lo_freq) +'MHz_noisy.fits'	
			data_quiet_name = 'cycle2_' + str(float(int(l))) +'_degrees_'+ str(lo_freq) +'MHz_quiet.fits'		
			
			print ('Galactic longtitude',l)			
			print (data_noisy_name)	
			print (data_quiet_name)
			
			data_on = Read_leushner_header(str(data_noisy_name))
			data_off = Read_leushner_header(str(data_quiet_name))
			H = Data_header(str(data_quiet_name))

			
			gain = 300 / np.sum(Average_spectrum(data_on, N=10, polarization_number='first') - Average_spectrum(data_off, N=10, polarization_number='first')) * np.sum(Average_spectrum(data_off, N=10, polarization_number='first'))
			
#			print (Average_spectrum(data_on, N=10, polarization_number='first'))
#			print (Average_spectrum(data_off, N=10, polarization_number='first'))

			
			ra, dec = H["RA"], H["DEC"]					
			print ('ra = ', int(ra), ' dec = ', int(dec))	

			N = H['NSPEC'] #number spectra
			print (N)
		
			
			doppler_velocity, Temperature  = final_calibration(data_quiet_name, Polarization_number='first', Gain= gain, velocity_correction='on' )
			
			#plot doppler vs temeprature: after longtitude 100 a weird graph appears
			plt.plot(doppler_velocity[2000:5000], Temperature[2000:5000] )
#			plt.xlim(-1000,1000)

#			print (len(doppler_velocity), len(Temperature))	
#			print ('doppler_velocity',max(doppler_velocity), min(doppler_velocity))
#			print ('Temperature',max(Temperature), min(Temperature))
			
						
			index_max = np.argmax(doppler_velocity[2000:5000])
			
			index = np.where(Temperature == max(Temperature[0:5500]))
			
#			print(index_max)
			print ('HI_temperature = ', max(Temperature), np.where(Temperature == max(Temperature)),  'HI_Velocity = ', doppler_velocity[index[0]][0])
#			peaks, _ = find_peaks(Temperature, height=0)	
#			print (peaks, max(Temperature[peaks]))
			#Plot average plots of all data in the file within the longtitude we collected
#			plot_average(Lo1_noisy_name, polarization_number= 'first')

			#Plot Temperatude vs Doppler velocity of all data in the file within the longtitude we collected
#			plot_temperature_doppler(final_data, velocity_correction= 'off')
			
			m = Basemap(projection='cyl') # Cylindrical projection
		
			Ra, Dec = m.makegrid(x_size , y_size)
			print('Ra', Ra)
			print('Dec', Dec)
			Ra_list = np.arange(0,360,1) # draw longitude lines every 30 degrees
			#lat_lines = np.arange(14,85,30) # draw latitude lines every 30 degrees
			dec_list = np.arange(-90,90,1)
			
#			Matrix = [[1 for x in np.arange(0,5,1)] for y in np.arange(-90,-80,1)] 
#			
#			Matrix[0][2]=23
			
#			plt.plot(doppler_velocity[2000:5000])
#			plt.plot(Temperature[2000:5000] )
#			plt.show()
			
			Grid[int(dec)][int(ra)]= max(Temperature[0:5500])
			Grid_velocity[int(dec)][int(ra)]= doppler_velocity[index[0]][0]
		

			plt.draw()
			plt.pause(.1)
			plt.close()
			
	plt.ioff()	
	
	plt.close()
	
	x = np.roll(Grid, 90, axis=0)
	x_vel = np.roll(Grid_velocity, 90, axis=0)
		
	plt.imshow(x, cmap='inferno', interpolation='gaussian', extent=[0,360,-90,90])
	
#	.drawmeridians(Ra_list, labels=[0,0,0,1], labelstyle='+/-')
#	m.drawparallels(dec_list, labels=[0,1,0,0], labelstyle='+/-')
	plt.colorbar()
	
	plt.xlabel('Right Ascension (degree)')
	plt.ylabel('Declination (degree)')
	plt.title("Temperature")

	plt.show()
	
	plt.imshow(x_vel, cmap='plasma', interpolation='gaussian', extent=[0,360,-90,90])
#	m.drawmeridians(Ra_list, labels=[0,0,0,1], labelstyle='+/-')
#	m.drawparallels(dec_list, labels=[0,1,0,0], labelstyle='+/-')
	plt.colorbar()
	plt.xlabel('Right Ascension (degree)')
	plt.ylabel('Declination (degree)')
	plt.legend()
	plt.show()

	y = np.roll(x, 180, axis=1)
	y_vel = np.roll(x_vel, 180, axis=1)
	
	plt.imshow(y, cmap='inferno', interpolation='gaussian', extent=[360,0,-90,90])
#	m.drawmeridians(Ra_list, labels=[0,0,0,1], labelstyle='+/-')
#	m.drawparallels(dec_list, labels=[0,1,0,0], labelstyle='+/-')
	plt.colorbar()
	plt.xlabel('Right Ascension (degree)')
	plt.ylabel('Declination (degree)')
	plt.title("Temperature")

	plt.show()
	
	plt.imshow(y_vel, cmap='inferno', interpolation='gaussian',extent=[360,0,-90,90])
#	m.drawmeridians(Ra_list, labels=[0,0,0,1], labelstyle='+/-')
#	m.drawparallels(dec_list, labels=[0,1,0,0], labelstyle='+/-')
	plt.colorbar()
	plt.xlabel('Right Ascension (degree)')
	plt.ylabel('Declination (degree)')
	plt.title("Velocity")

	plt.show()
	
	





########################## Testing old data ##############################################################

plot_all_attempts()
fon= fits.open('attempt_all_160.0_degrees_635MHz_noisy.fits')
N = fon[0].header['NSPEC']
	

foff= fits.open('attempt_all_160.0_degrees_635MHz_noisy.fits')
N = foff[0].header['NSPEC']

print (foff[0].header)
gain = Gain(fon, foff, N, 301, 1)
#calculate the gain:
print ('Gain', Gain(fon, foff, N, 301, 1))


#comparing on/off data, by inputing file name and polarization number 

#plot_average('twenty_plane_on.fits')
#plot_average('twenty_plane_off.fits')
#plt.show()
#
#plot_average('twenty_plane_on.fits','second')
#plot_average('twenty_plane_off.fits','second')
#plt.show()



#comparing on/off data, by inputing file name and polarization number 

#plot_average('attempt_all_160.0_degrees_635MHz_noisy.fits')
#plot_average('attempt_all_160.0_degrees_634MHz_noisy.fits')

#plot_average('attempt_all_160.0_degrees_635MHz_quiet.fits')
#plot_average('attempt_all_160.0_degrees_634MHz_quiet.fits')
#
#plot_average('attempt_all_160.0_degrees_635MHz_quiet.fits', 'second')
#plot_average('attempt_all_160.0_degrees_634MHz_quiet.fits', 'second')


#plot_average('attempt_all_160.0_degrees_635MHz_noisy.fits')
#plot_average('attempt_all_160.0_degrees_634MHz_noisy.fits')

#plot_average('attempt_all_160.0_degrees_635MHz_noisy.fits', 'second')
#plot_average('attempt_all_160.0_degrees_634MHz_noisy.fits', 'second')
#plt.show()

#plot_average('attempt_all_160.0_degrees_634MHz_noisy.fits','second')
#plot_average('attempt_all_160.0_degrees_635MHz_noisy.fits','second')

#plt.show()

#plot all rawa data from file
#Plot_leushner_raw_data('thirty_plane_on.fits')









def plot_temperature_doppler(data, velocity_correction= 'off'):
	
	
	data_x = data[0]
	data_y = data[1]			
	
	plt.plot(data_x, data_y, label=' Doppler correction is '+ velocity_correction)
	
	plt.xlabel('Doppler velocity (km/s)')
	plt.ylabel('Temperature (kelvin)')
	plt.legend()
	


#calibration, (ugradio.doppler.get_projected_velocity(ra, dec, jd)) outputs the doppler velocity and the unit, this format blocks the velocity correction calculation when it's on.


#data = final_calibration('attempt_all_160.0_degrees_635MHz_quiet.fits', velocity_correction='off', Gain= 2000)
#plt.plot(data[0], data[1], label='No Doppler correction')
#plt.xlabel('Doppler velocity (km/s)')
#plt.ylabel('Temperature (kelvin)')
#plt.legend()
#plt.show()	
			
#plot_all_attempts(634, 635, file_name = 'attempt_all')	
#plt.show()	

	



#data = final_calibration('attempt_all_160.0_degrees_635MHz_quiet.fits', velocity_correction='on', Gain= gain)
#plt.plot(data[0], data[1], label=' Doppler corrected' )
#plt.xlabel('Doppler velocity (km/s)')
#plt.ylabel('Temperature (kelvin)')
#plt.legend()
#plt.show()
#
