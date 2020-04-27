
import matplotlib.pyplot as plt
import numpy as np
import scipy
import time

import ugradio.leusch
import ugradio.timing
import ugradio.leo as leo

from astropy import units as u
from astropy.time import Time
import astropy.io.fits as fits
from astropy.coordinates import SkyCoord,AltAz,EarthLocation


#Altitude boundaries in degrees
minimum_alt, maximum_alt = 15, 85

#Azimuth boundaries in degrees
minimum_az, maximum_az = 5, 350

#Parameters 
degree_spacing = 2
lontitude_min= -10
longtitude_max = 250

Longtitude_Spacing = (longtitude_max - (lontitude_min)) / degree_spacing
	

#setting up the telescope coordinates on earth from astropy earthlocation
telescope_coordinates = EarthLocation(lat=ugradio.leo.lat*u.deg, lon=ugradio.leo.lon*u.deg, height=ugradio.leo.alt*u.m) 

print ("telescope_coordinates (lat = {0} , lon = {1}, height= {2} ) ".format(telescope_coordinates.value[0],telescope_coordinates.value[1],telescope_coordinates.value[2]))


def longitudes_values_Missing(file_name='attempt_all'):
	
	'''
	find the missing longtitudes from the data files we already took using the file name,
	
	the assumption is that if we don't find the longtitude value in the file name means that we didn't take that data at that specific longtitude 
	'''
	
	degree_spacing = 2/2

	lontitude_min= -10
	longtitude_max = 250
	Longtitude_Spacing = (longtitude_max - (lontitude_min)) / degree_spacing
	
	longitude_range = np.linspace(lontitude_min, longtitude_max, int(Longtitude_Spacing)+1)	
	
	Missing_longitudes_values = []
	
	for l in longitude_range:
		
		try:				
#			Check if langitude is taken already or not	
			fits.open(file_name + '_' + str(float(int(l))) + '_degrees_635MHz_noisy.fits')
		except:
			Missing_longitudes_values.append(l)
			
	np.save("missing_longitudes_values", Missing_longitudes_values)
	
	return Missing_longitudes_values
	

#take on and off noise data that we will use for temperature calibration at the end 
def Capture_noise_on_and_off(Lo_freq, N, l, b,  file_name = 'attempt_all'):
	
	spectrometer = leuschner.Spectrometer('10.0.1.2')
	noise = ugradio.leusch.LeuschNoise()
	
	noise.on()
	spectrometer.read_spec(file_name + '_' + str(l)+ '_degrees_' + Lo_freq + 'MHz_noisy.fits', N, (l, b), 'ga')
	
	noise.off()
	spectrometer.read_spec(file_name + '_' + str(l)+ '_degrees_' + Lo_freq + 'MHz_quiet.fits', N, (l, b), 'ga')
	

	
def Galactic_to_Topocentric_converter(l, b, jd= ugradio.timing.julian_date()):
	
	'''
	Converts galactic coordinates longitude l and latitude b to topocentric coordinates altitude and azimuth 		
	We will use the altitude and azimuth to poin the telescope toward the galacric plane where b=0 and -10 <l< 250
	'''
	

	#active when we are taking data from telecope
	#not active when we want to convert l/b to alt/az in a different julian date 
#	jd= ugradio.timing.julian_date()
	
	#set up the galatic coordinates	using astropy skycood frame
	Galactic_position = SkyCoord(frame='galactic', l=l, b=b, unit=(u.degree,u.degree))
	
	#transform into altitude an azimuth using astropy AltAz in degrees
	Alt_az= Galactic_position.transform_to(AltAz(obstime=Time(jd,format='jd'),location=telescope_coordinates))
	
	altitude, azimuth = Alt_az.alt.deg, Alt_az.az.deg	
	
	print('altitude, azimuth',altitude, azimuth)
	
	return altitude, azimuth


#print (Galactic_to_Topocentric_converter(200,0)[0])



def Capture_galactic_data(file_name, Lo_freq1, Lo_freq2, N=10):
	
	'''	
	Capture long term data using the parameters longtidude range mentioned bellow
		
	'''
	
	#parameters
	degree_spacing = 2/2

	lontitude_min= -10
	longtitude_max = 250
	Longtitude_Spacing = (longtitude_max - (lontitude_min)) / degree_spacing
	
	
	b = 0 #galactic plane	
	N = 10 # number of spectra taking at each longtitude, we shold probably increase it to 15 or 20, it would make the average better. 
	longitude_range = np.linspace(lontitude_min, longtitude_max, int(Longtitude_Spacing)+1)	 # range of longtitude
	
	#initiating leushner telecope
	LO = ugradio.agilent.SynthDirect()
	LT = ugradio.leusch.LeuschTelescope()
	spectrometer = leuschner.Spectrometer('10.0.1.2')
	
	
	
	

#get the data of different longtitudes that were not taking before to complete a larger map
	for l in longitudes:	
								
		alt, az = Galactic_to_Topocentric_converter(l, b) # converts l/b to alt/az
		
		#makes sure the alt, az is within the telescopes boundries
		if (minimum_alt < alt < maximum_alt) and (minimum_az < az < maximum_az): 
			
			LT.point(alt, az)
			
			#set LO frequency
			LO.set_frequency(Lo_freq1, "MHz")
			#take on and off noise data that we will use for temperature calibration 
			Capture_noise_on_and_off( Lo_freq1, N, l, b,file_name = 'attempt_all' )
			
			
			#set LO frequency
			LO.set_frequency(Lo_freq2, "MHz")
			#take on and off noise data that we will use for temperature calibration 
			Capture_noise_on_and_off( Lo_freq2, N, l, b, file_name = 'attempt_all' )
	

	
	
		
	LT.stow()
	
	#save the values of the longtitude were data is missing
	longitudes_values_Missing(file_name)


def Capture_missing_galactic_data(file_name, Lo_freq1, Lo_freq2, missing_data_file_name = "missing_longitudes_values", N=10):
	
	'''	
	Capture long term missing data using the parameters longtidude range mentioned bellow
	'''
	
	#parameters
	degree_spacing = 2/2

	lontitude_min= -10
	longtitude_max = 250
	Longtitude_Spacing = (longtitude_max - (lontitude_min)) / degree_spacing
		
		
	b = 0 #galactic plane	
	N = 10 # number of spectra taking at each longtitude, we shold probably increase it to 15 or 20, it would make the average better. 
	longitude_range = np.linspace(lontitude_min, longtitude_max, int(Longtitude_Spacing)+1)	 # range of longtitude
	
	
	LO = ugradio.agilent.SynthDirect()
	LT = ugradio.leusch.LeuschTelescope()
	spectrometer = leuschner.Spectrometer('10.0.1.2')
	
	#set LO frequency
	LO.set_frequency(635, "MHz")
	
	#get the corresponding longtitude of the missing data sets	
	missing_longitudes = np.load(missing_data_file_name)


#get the data of different longtitudes that were not taking before to complete a larger map
	for l in longitudes:
		
		# makes sure the data from the corresponding longitude is missing to take new data
		if l in missing_longitudes:
						
			alt, az = Galactic_to_Topocentric_converter(l, b) # converts l/b to alt/az
			
			#makes sure the alt, az is within the telescopes boundries
			if (minimum_alt < alt < maximum_alt) and (minimum_az < az < maximum_az): 
				
				LT.point(alt, az)
				
				#set LO frequency
				LO.set_frequency(Lo_freq1, "MHz")
				#take on and off noise data that we will use for temperature calibration 
				Capture_noise_on_and_off(Lo_freq1, N, l, b, file_name = 'attempt_all')
				
				
				#set LO frequency
				LO.set_frequency(Lo_freq2, "MHz")
				#take on and off noise data that we will use for temperature calibration 
				Capture_noise_on_and_off(Lo_freq2 , N, l, b , file_name = 'attempt_all' )
				

	
	#take on and off noise data that we will use for temperature calibration at the end 
	
		
	LT.stow()
	
	#save the values of the longtitude were data is missing
	longitudes_values_Missing(file_name)


	
