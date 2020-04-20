
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
Longtitude_Spacing = (250 - (-1*10)) / degree_spacing 

#setting up the telescope coordinates on earth from astropy earthlocation
telescope_coordinates = EarthLocation(lat=ugradio.leo.lat*u.deg, lon=ugradio.leo.lon*u.deg, height=ugradio.leo.alt*u.m) 

print ("telescope_coordinates (lat = {0} , lon = {1}, height= {2} ) ".format(telescope_coordinates.value[0],telescope_coordinates.value[1],telescope_coordinates.value[2]))


def longitudes_values_Missing(file_name):
	
	'''
	find the missing longtitudes from the data files we already took using the file name,
	
	the assumption is that if we don't find the longtitude value in the file name means that we didn't take that data at that specific longtitude 
	'''
	
	degree_spacing = 2

	lontitude_min= -10
	longtitude_max = 250
	Longtitude_Spacing = (250 - (-1*10)) / degree_spacing
	
	longitude_range = np.linspace(lontitude_min, longtitude_max, Longtitude_Spacing)
	
	Missing_longitudes_values = []
	
	for l in longitude_range:
		
		try:		
				
			fits.open("Milkyway/" + file_name + str(l) + ".fits")
			
		except:
			Missing_longitudes_values.append(l)
			
	np.save("missing_longitudes_values", Missing_longitudes_values)
	
	return Missing_longitudes_values
	

#take on and off noise data that we will use for temperature calibration at the end 
def Capture_noise_calibration(file_name, N, l, b):
	
	spectrometer = leuschner.Spectrometer('10.0.1.2')
	noise = ugradio.leusch.LeuschNoise()
	
	noise.on()
	spectrometer.read_spec("Leushner/" + file_name + str(l) + "_noise_on.fits", N, (l, b), 'ga')
	
	noise.off()
	spectrometer.read_spec("Leushner/" + file_name + str(l) + "_noise_off.fits", N, (l, b) ,'ga')
	

	
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
	
	print(altitude, azimuth)
	
	return altitude, azimuth


print (Galactic_to_Topocentric_converter(200,0)[0])

def Capture_galactic_data(file_name, N=10):
	
	'''	
	Capture long term data using the parameters longtidude range mentioned bellow
		
	'''
	
	#parameters
	degree_spacing = 2
	lontitude_min= -10
	longtitude_max = 250
	Longtitude_Spacing = (250 - (-1*10)) / degree_spacing		
	b = 0 #galactic plane	
	N = 10 # number of spectra taking at each longtitude, we shold probably increase it to 15 or 20, it would make the average better. 
	longitudes = np.linspace(lontitude_min, longtitude_max, Longtitude_Spacing) # range of longtitude
	
	#initiating leushner telecope
	LO = ugradio.agilent.SynthDirect()
	LT = ugradio.leusch.LeuschTelescope()
	spectrometer = leuschner.Spectrometer('10.0.1.2')
	
	
	#set LO frequency
	LO.set_frequency(635, "MHz")
	


#get the data of different longtitudes that were not taking before to complete a larger map
	for l in longitudes:	
								
		alt, az = Galactic_to_Topocentric_converter(l, b) # converts l/b to alt/az
		
		#makes sure the alt, az is within the telescopes boundries
		if (minimum_alt < alt < maximum_alt) and (minimum_az < az < maximum_az): 
			
			LT.point(alt, az)
			
			#save data in a file called Leushner  
			spectrometer.read_spec("Leushner/" + file_name + str(l) + ".fits", N, (l, b), 'ga')
				

	
	#take on and off noise data that we will use for temperature calibration at the end 
	Capture_noise_calibration(file_name, N, l, b)
		
	LT.stow()
	
	#save the values of the longtitude were data is missing
	longitudes_values_Missing(file_name)


def Capture_missing_galactic_data(file_name, missing_data_file_name = "missing_longitudes_values", N=10):
	
	'''	
	Capture long term missing data using the parameters longtidude range mentioned bellow
	'''
	
	#parameters
	degree_spacing = 2
	lontitude_min= -10
	longtitude_max = 250
	Longtitude_Spacing = (250 - (-1*10)) / degree_spacing		
	b = 0 #galactic plane	
	N = 10 # number of spectra taking at each longtitude, we shold probably increase it to 15 or 20, it would make the average better. 
	longitudes = np.linspace(lontitude_min, longtitude_max, Longtitude_Spacing) # range of longtitude
	
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
				
				#save data in a file called Leushner  
				spectrometer.read_spec("Leushner/" + file_name + str(l) + ".fits", N, (l, b), 'ga')
				

	
	#take on and off noise data that we will use for temperature calibration at the end 
	Capture_noise_calibration(file_name, N, l, b)
		
	LT.stow()
	
	#save the values of the longtitude were data is missing
	longitudes_values_Missing(file_name)


	
