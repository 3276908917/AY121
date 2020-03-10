import numpy as np
import scipy as sp
import ugradio
import matplotlib.pyplot as plt
import time

class Interferometery:
	
	def __init__(self, Update_position_time , julian_day, latitude, longitude, altitude ):
		self.Update_position_time = Update_position_time
		self.julian_day = julian_day
		self.latitude = latitude
		self.longitude = longitude
		self.altitude = altitude
		
		
	def initialize_control(self):	
		return ugradio.interf.interferometer()
		
	def initialize_voltage(self):
		return ugradio.hp_multi.HP_Multimeter()

	def power_spectrum(self,signal):
		fourier_transform = np.fft.fft(signal)
		return np.abs(fourier_transform)**2


	def fourier_frequency(self, signal, sampling_frequency):
		
		signal_fft=np.fft.fft(signal)

		time_step = 1 / sampling_frequency 
		freqs = np.fft.fftfreq(signal_fft.size, time_step)
		index = np.argsort(freqs)
		
		plt.plot(freqs[index]/1e6, ps[index])
		plt.xlabel('Frequency (MHz)')
		plt.ylabel('Amplitude (V)')
		plt.title(' Fourier transform ')
#		plt.yscale('log')
		
		plt.legend()
		plt.show()


		return freqs[index], signal_fft[index]
		
	def spectrum_frequency(self, signal, sampling_frequency):

		ps=power_spectrum(signal)
		time_step = 1 / sampling_frequency 
		freqs = np.fft.fftfreq(signal.size, time_step)
		index = np.argsort(freqs)
		
		#plot the frequency power spectrum
		plt.plot(freqs[index]/1e6, ps[index])
		plt.xlabel('Frequency (MHz)')
		plt.ylabel('Amplitude (V^2)')
		plt.title(' Power spectrum ')
#		plt.yscale('log')
		
		plt.legend()
		plt.show()

		return freqs[index], ps[index]
	
	def Test_system(self, sampling_frequency):
		
		ifm = initialize_control(self)
		hpm = initialize_voltage(self)
		
		ifm.maintenance()
		
		
		data = hpm.read_voltage()
		
		Interferometery.fourier_frequency(self, data, sampling_frequency)
						
		return np.save('test_data', data)
		
		
	def Record_sun(self, recording_time, time_delay=0):
		
		
		ifm = initialize_control(self)
		hpm = initialize_voltage(self)
		
		time_1 = time.time()
		
		while self.Update_position_time >= time.time() - time_1:		
			
			count = 0
			
			Ra, Dec = ugradio.coord.sunpos([self.julian_day])
			
			Alt, azimuth = ugradio.coord.get_altaz(Ra, Dec [ self.julian_day [self.latitude [ self.longitude [ self.altitude [ equinox ]]]]])
			
			ifm.point(Alt, azimuth)
			
			initial_time = time.time()
					
			
			while recording_time >= time.time() - initial_time :
				
				time.sleep(time_delay)				
				hpm.start_recording(recording_time)
			
									
			data = hpm.get_recording_data()
			
			data_name = 'data_sun_'+ str(count)
					
			np.save( data_name, data)		
			
					
			hpm.end_recording()
			
			
			count+=1
		
			
		ifm.stow()
			
		return count
			
	def Record_moon(self, recording_time, time_delay=0):
		

		ifm = initialize_control(self)
		hpm = initialize_voltage(self)
		
		time_1 = time.time()
		
		while time_1 - time.time() <= self.Update_position_time:		
			
			count = 0
			
			Ra, Dec = ugradio.coord.moonpos([ self.julian_day [ self.latitude [ self.longitude [ self.altitude  ]]]])
			
			Alt, azimuth = ugradio.coord.get_altaz(Ra, Dec [ self.julian_day [self.latitude [ self.longitude [ self.altitude [ equinox ]]]]])
			
			ifm.point(Alt, azimuth)
			
			initial_time = time.time()
					
			
			
			while time.time() <= initial_time + recording_time:	
				
				time.sleep(time_delay)				
				hpm.start_recording(recording_time)
							
			
			data = hpm.get_recording_data()
			
			data_name = 'data_moon_'+ str(count)
					
			np.save( data_name, data)		
			
					
			hpm.end_recording()
			
			
			count+=1
		
			
		ifm.stow()
			
		return count
		
	def Record_star(self, ra, dec, equinox, recording_time, total_observation_time):
		
		count=0
		
		ifm = initialize_control(self)
		hpm = initialize_voltage(self)
		
		Ra, Dec = ugradio.coord.precess(ra, dec [ self.julian_day [ equinox ]])
		
		Alt, azimuth = Alt, azimuth = ugradio.coord.get_altaz(Ra, Dec [ self.julian_day [self.latitude [ self.longitude [ self.altitude [ equinox ]]]]])
		
		ifm.point(Alt, azimuth)

		intitial_time = time.time()		
		
		while   time.time() <= total_observation_time + intitial_time:
			
			print(ifm.get_pointing())
				
			while time.time() <= initial_time + recording_time:
				
				time.sleep(time_delay)				
				hpm.start_recording(recording_time)
											
				
			data = hpm.get_recording_data()
			
			data_name = 'data_moon_'+ str(count)
					
			np.save( data_name, data)		
			
					
			hpm.end_recording()
			
			
			count+=1
		
			
		
		ifm.stow()
			
		return count
				
			
