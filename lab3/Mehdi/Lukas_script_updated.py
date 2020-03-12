import ugradio
import time
import numpy as np

## aliasing

def jd():
	return ugradio.timing.julian_date()

class Irf:
	def __init__(self, equinox='J2000',
		latitude=ugradio.coord.nch.lat, longitude=ugradio.coord.nch.lon, altitude=ugradio.coord.nch.alt,
	):
		self.lat = latitude
		self.lon = longitude
		self.alt = altitude
		self.eq = equinox

		self.ctrl = ugradio.interf.Interferometer()
		self.multi = ugradio.hp_multi.HP_Multimeter()

	def test_system(self):
		self.ctrl.maintenance()

		data = self.multi.read_voltage()
		print(data)

		self.ctrl.stow()
		return np.savez('test_data', data)

	# Idea: query user for inputs, then run the data collection routine

	### Section: positioning functions. Get altitude and azimuth for...

	def sun(self):
		ra_sun, dec_sun = ugradio.coord.sunpos(jd())
		return ugradio.coord.get_altaz(ra_sun, dec_sun, jd(), self.lat, self.lon, self.alt, self.eq)
	
	def moon(self):
		ra_moon, dec_moon = ugradio.coord.moonpos(jd(), self.lat, self.lon, self.alt)
		return ugradio.coord.get_altaz(ra_moon, dec_moon, jd(), self.lat, self.lon, self.alt, self.eq)

	def star(self, old_ra, old_dec):
		
		'''
		I am not sure about my syntax here, it looks a little fishy.
		Specifically, I am not sure if the child function will keep track of the variables.
		
		
		- stargazer will keep track of all the variables and will still work perfectly. But it's useless to me, unless you had a reason for it??
		
		'''
		
		
#		def stargazer():
#			new_ra, new_dec = ugradio.coord.precess(old_ra, old_dec, jd(), self.eq)
#			return ugradio.coord.get_altaz(new_ra, new_dec, jd(), self.lat, self.lon, self.alt, self.eq)
#		return stargazer()
			
		new_ra, new_dec = ugradio.coord.precess(old_ra, old_dec, jd(), self.eq)
		return ugradio.coord.get_altaz(new_ra, new_dec, jd(), self.lat, self.lon, self.alt, self.eq)
			
		
		
	
	def coord(self, name='sun', old_ra=50, old_dec=50):	
			
		if name=='sun' or name=='Sun' :
			return irf.sun(self)
			
		elif name=='moon' or name=='Moon':
			return irf.moon(self)
			
		elif name=='star' or name=='Star':
			return irf.star(self, old_ra, old_dec)
			
			


	### end section
	def point(self, name='sun', old_ra=50, old_dec=50)):
		
		alt_target, az_target = irf.coord(self, name, old_ra=50, old_dec=50))
		
		self.ctrl.point(alt_target, az_target, wait = True)
		
		actual = self.ctrl.get_pointing()

		if  175 > actual['ant_w'][0] > 5 \
			or 175 > actual['ant_e'][0] > 5  \
			or 300 > actual['ant_w'][1] > 90  \
			or 300 > actual['ant_e'][1] > 90:					
			raise AssertionError('Target is out of range!')
		
		print ( 'Actual coordinates' ,	actual )
		return True
			



	def capture(self, label,
				total_capture_time = 3960, reposition_interval = 60,
				backup_interval = 600, capture_interval = 1):
		'''
		The telecopes will start to collect data before pointing at the source, it's ok, we can just eliminate the first minute of data collection
		'''
		
		recording_start = last_backup = pointng_time = time.time()

		index = 0
		
		irf.point(self)
		
		self.multi.start_recording(capture_interval)

		while total_capture_time >= time.time() - recording_start :
			
			
			"""
				
			Remarks:
				
			- Using your setup , the telecopes will point continuously at the target and only stops pointing when we collect data every 10 minutes. 
			
			- time.sleep(reposition_interval) will stop the code for 60 seconds. If time.sleep doesn't interupt collecting data then it's ok. if not, it will be collecting only one data set after it's pointed. Also just noticed that get_pointing had wait=true, so 
			
			- Using backup_interval and last_backup and resetting the last_backup to the actual time was a good idea. 
			
			- I wonder if it will save accumulated data instead if we don't end the recording hpm.end_recoding(). But in your setup, if we end the recording it won't save anything anymore, since you called it before the while loop.   ????
			
		
			"""
			
			
			
			# if time.sleep() interupts with the data collection try this instead, if not you can keep time.sleep(reposition_interval)
			
			if time.time() - pointng_time >= reposition_interval:
				
				irf.point(self)
				
				pointng_time=time.time()
				
				time.sleep(3)
								
#			irf.point(self)
				
			if time.time() - last_backup >= backup_interval:
				data = self.multi.get_recording_data()
				# create a time-stamp for the data
				minutes = str(np.around((time.time() - recording_start) / 60, 2))
				data_name = 'data/' + label + '_' + minutes + '_minutes'
				np.savez(data_name, data=data)
				last_backup = time.time()
				
				self.multi.end_recording()
				self.multi.start_recording(capture_interval)
				
#			time.sleep(reposition_interval)
	
	
	
	
			

		self.ctrl.stow()
		print('No runtime errors encountered.')