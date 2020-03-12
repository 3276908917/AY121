import numpy as np
import scipy as sp
import ugradio
import matplotlib.pyplot as plt
import time

class Interferometery:
	
	def __init__(
                self, Total_recording_time, Update_position_time, data_saved_time,
                julian_day=ugradio.timing.julian_date(), latitude=ugradio.coord.nch.lat,
                longitude=ugradio.coord.nch.lon, altitude=ugradio.coord.nch.alt, equinox='J2000'):
                # input all the location informations
	
		self.Total_recording_time= Total_recording_time  # total reading time, 1.1 hour in our case
		self.Update_position_time = Update_position_time  # time to switch the pointing position of the telecopes, 1 min is good i guess
		self.julian_day = julian_day 
		self.latitude = latitude
		self.longitude = longitude
		self.altitude = altitude
		self.equinox = equinox
		self.data_saved_time = data_saved_time # time to save every data set, 10 mins is good i guess
		
	def initialize_control(self):	
		return ugradio.interf.Interferometer()
		
	def initialize_voltage(self):
		return ugradio.hp_multi.HP_Multimeter()

	def power_spectrum(self, signal):
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
		
		ifm = self.initialize_control()
		hpm = self.initialize_voltage()
		
		ifm.maintenance()
		
		
		data = hpm.read_voltage()
		
		print (data)
#		Interferometery.fourier_frequency(self, data, sampling_frequency)
						
		return np.savez('test_data', data)
		
		
	def Record_sun(self, recording_time):
		
		
		ifm = self.initialize_control()
		hpm = self.initialize_voltage()
		initial_time = time.time()
		
		index = 0
		hpm.start_recording(recording_time)
		
		while self.Total_recording_time >= time.time() - initial_time : # it will read for an hour if total recording time is an hour
			
			Ra, Dec = ugradio.coord.sunpos(self.julian_day) # get's RA and dec of the sun
			
			Alt, azimuth = ugradio.coord.get_altaz(Ra, Dec, self.julian_day, self.latitude, self.longitude, self.altitude, self.equinox)
			
			ifm.point(Alt, azimuth, wait=True) # we can add (wait=True) to wait until it's pointed to proceed
			
			#time.sleep(2) # to give our telecopes time to point at the sun


			if self.data_saved_time * index <= time.time() - initial_time:
					
				data = hpm.get_recording_data()
					
				data_name = 'data/data_sun_' + str(index + 1)
							
				np.savez(data_name, data)		
							
				hpm.end_recording()
					
				index += 1

			time_1 = time.time()
							
			while self.Update_position_time >= time.time() - time_1: # it will read data until it's time to switch position ( every 1 minute is better i guess) 
				time.sleep(1)                        								
				 # record 1 data every 'recording_time' seconds
				
							
#			the data will be saved every time the telescope change position (1 minute) not sure if it's convinient, probably 10 mins is better, so i concidered saving the data every 10th time that we change our position. 
			
		ifm.stow()
			
		return index
	
	def Record_moon(self, recording_time):
		
		
		ifm = self.initialize_control()
		hpm = self.initialize_voltage()
		initial_time = time.time()
		
		a=1
		
		while self.Total_recording_time >= time.time() - initial_time : # it will read for an hour if total recording time is an hour
		
			count = 0
			
			Ra, Dec = ugradio.coord.moonpos(self.julian_day, self.latitude, self.longitude, self.altitude)
			
			Alt, azimuth = ugradio.coord.get_altaz(Ra, Dec, self.julian_day, self.latitude, self.longitude, self.altitude, self.equinox)
			
			ifm.point(Alt, azimuth, wait=True) # we can add (wait=true) to wait until it's pointed to proceed
			
			#time.sleep(2) # to give it time to point 
			

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
		
					
	def Record_star(self, old_ra, old_dec, recording_time, total_observation_time):
		
		count = 0
		
		ifm = self.initialize_control()
		hpm = self.initialize_voltage()
		
		new_ra, new_dec = ugradio.coord.precess(old_ra, old_dec, self.julian_day, self.equinox)
		
		Alt, azimuth = ugradio.coord.get_altaz(new_ra, new_dec, self.julian_day, self.latitude, self.longitude, self.altitude, self.equinox)
		
		ifm.point(Alt, azimuth)

		intitial_time = time.time()		
		
		while time.time() <= total_observation_time + intitial_time:
			
			print(ifm.get_pointing())
				
			while time.time() <= initial_time + recording_time:
				
				time.sleep(time_delay)				
				hpm.start_recording(recording_time)
												
			data = hpm.get_recording_data()
			
			data_name = 'data_moon_'+ str(count)
					
			np.save(data_name, data)		
			
			hpm.end_recording()
			
			count += 1
		
		ifm.stow()
			
		return count
				
			
