import numpy as np
import scipy as sp
import ugradio
import matplotlib.pyplot as plt
import time

class Interferometery:
	
	def __init__(self, Total_recording_time, Update_position_time , data_saved_time, julian_day, latitude, longitude, altitude, equinox ): # input all the location informations
	
		self.Total_recording_time= Total_recording_time
		self.Update_position_time = Update_position_time
		self.julian_day = julian_day
		self.latitude = latitude
		self.longitude = longitude
		self.altitude = altitude
		self.equinox = equinox
		self.data_saved_time = data_saved_time
		
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
		
		print (data)
#		Interferometery.fourier_frequency(self, data, sampling_frequency)
						
		return np.savez('test_data', data)
		
		
	def Record_sun(self, recording_time, time_delay=0):
		
		
		ifm = initialize_control(self)
		hpm = initialize_voltage(self)
		initial_time = time.time()
		
		a=1
		
		while self.Total_recording_time >= time.time() - initial_time : # it will read for an hour if total recording time is an hour
		
			count = 0
			
			Ra, Dec = ugradio.coord.sunpos([self.julian_day])
			
			Alt, azimuth = ugradio.coord.get_altaz(Ra, Dec [ self.julian_day [self.latitude [ self.longitude [ self.altitude [ self.equinox ]]]]])
			
			ifm.point(Alt, azimuth) # we can add (wait=true) to wait until it's pointed to proceed
			
			time.sleep(2) # to give it time to point 
			

			time_1 = time.time()
			
							
			while self.Update_position_time >= time.time() - time_1: # it will read data untill it's time to switch position ( every 1 minute is better ) 
											
				hpm.start_recording(recording_time) 
				
							
#			the data will be saved every time the telescope change position (1 minute) not sure if it's convinient, probably 10 mins is better, so i concidered saving the data every 10th time that we change our position. 
			
				if 	self.data_saved_time * a >= time.time() - initial_time:
					
					data = hpm.get_recording_data()
					
					data_name = 'data_sun_' + str(count)
							
					np.savez(data_name, data)		
							
					hpm.end_recording()
					
					a+=1
			
		
					count+=1
		
			
		ifm.stow()
			
		return count
	
	def Record_moon(self, recording_time, time_delay=0):
		
		
		ifm = initialize_control(self)
		hpm = initialize_voltage(self)
		initial_time = time.time()
		
		a=1
		
		while self.Total_recording_time >= time.time() - initial_time : # it will read for an hour if total recording time is an hour
		
			count = 0
			
			Ra, Dec = ugradio.coord.moonpos([ self.julian_day [ self.latitude [ self.longitude [ self.altitude  ]]]])
			
			Alt, azimuth = ugradio.coord.get_altaz(Ra, Dec [ self.julian_day [self.latitude [ self.longitude [ self.altitude [ self.equinox ]]]]])
			
			ifm.point(Alt, azimuth) # we can add (wait=true) to wait until it's pointed to proceed
			
			time.sleep(2) # to give it time to point 
			

			time_1 = time.time()
			
							
			while self.Update_position_time >= time.time() - time_1: # it will read data untill it's time to switch position ( every 1 minute is better ) 
											
				hpm.start_recording(recording_time) 
				
							
#			the data will be saved every time the telescope change position (1 minute) not sure if it's convinient, probably 10 mins is better, so i concidered saving the data every 10th time that we change our position. 
			
				if 	self.data_saved_time * a >= time.time() - initial_time:
					
					data = hpm.get_recording_data()
					
					data_name = 'data_sun_' + str(count)
							
					np.savez(data_name, data)		
							
					hpm.end_recording()
					
					a+=1
			
		
					count+=1
		
			
		ifm.stow()
			
		return count
			
			
					
	def Record_star(self, ra, dec, recording_time, total_observation_time):
		
		count=0
		
		ifm = initialize_control(self)
		hpm = initialize_voltage(self)
		
		Ra, Dec = ugradio.coord.precess(ra, dec [ self.julian_day [ self.equinox ]])
		
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
				
			
