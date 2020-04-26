import matplotlib.pyplot as plt
import numpy as np
import scipy
import time

import ugradio.leusch as leusch
import ugradio.timing
import ugradio.leo as leo

from astropy import units as u
from astropy.time import Time
import astropy.io.fits as fits
from astropy.coordinates import SkyCoord, AltAz, EarthLocation
 
degree_spacing = 2
# galactic coordinates and ranges
g_lat = 0
g_lon_min = -10
g_lon_max = 250
longtitude_spacing = (250 - (-1*10)) / degree_spacing

longitudes = np.linspace(
    g_lon_min,
    g_lon_max,
    longtitude_spacing
)

# L: what is this argument??
spectrometer = leuschner.Spectrometer('10.0.1.2')
LO = ugradio.agilent.SynthDirect()
LT = ugradio.leusch.LeuschTelescope()
#set LO frequency
LO.set_frequency(635, "MHz")

#setting up the telescope coordinates on earth from astropy earthlocation
telescope_coordinates = EarthLocation(
    lat=ugradio.leo.lat * u.deg,
    lon=ugradio.leo.lon * u.deg,
    height=ugradio.leo.alt * u.m
) 

print ("telescope_coordinates (lat = {0} , lon = {1}, height= {2} ) ".format(
    telescope_coordinates.value[0],
    telescope_coordinates.value[1],
    telescope_coordinates.value[2]
))

def longitudes_values_missing(file_name):
    '''
    Find the missing longtitudes from the data files we already took using the file name.
    The assumption is that if we don't find the longtitude value in the file name,
    we didn't take that data at that specific longtitude 
    '''
    missing_longitudes_values = []

    for l in longitudes:
        try:		                        
            fits.open("Plane/" + file_name + str(l) + ".fits")    
        except:
            missing_longitudes_values.append(l)
                    
    np.save("missing_longitudes_values", missing_longitudes_values)
    return missing_longitudes_values
	 
def capture_noise_calibration(file_name, N, l, b):
    '''
    take on and off noise data that we will use
    for temperature calibration at the end
    '''
    noise = ugradio.leusch.LeuschNoise()
    
    noise.on()
    spectrometer.read_spec(file_name + str(l) + "_noise_on.fits", N, (l, b), 'ga')
    
    noise.off()
    spectrometer.read_spec(file_name + str(l) + "_noise_off.fits", N, (l, b),'ga')

# Labeling conflict with my counterpart function!
def gal_to_topo(l, b, jd):
    '''
    Converts galactic coordinates longitude l and latitude b to
    topocentric coordinates altitude and azimuth.
    We use the alt and az to point the telescope
    toward the galactic plane where b = 0 and -10 < l < 250.
    '''
    #active when we are taking data from telecope
    #not active when we want to convert l/b to alt/az in a different julian date 
    jd= ugradio.timing.julian_date()
	
    #set up the galatic coordinates using astropy skycood frame
    galactic_position = SkyCoord(
        frame='galactic', l=l, b=b, unit=(u.degree,u.degree)
    )
    
    #transform into altitude an azimuth using astropy AltAz in degrees
    alt_az= galactic_position.transform_to(
        AltAz(obstime=Time(jd, format='jd'),
        location=telescope_coordinates)
    )
    
    alt, az = alt_az.alt.deg, alt_az.az.deg	
    print(alt, az)
    return alt, az

print (galactic_to_topocentric_converter(200, 0)[0])

# Mehdi suggests increasing N to 15 or 20
# to improve averaging.
def capture_plane_data(file_name, N=10):
    '''	
    Capture long term data using the
    universal longtidude range parameter
    (top of script)
    '''			
    #initiating leushner telecope
    LO = ugradio.agilent.SynthDirect()
    LT = ugradio.leusch.LeuschTelescope()
    #set LO frequency
    LO.set_frequency(635, "MHz")

    # get the data of different longtitudes that
    # we're not taking before to complete a larger map
    for l in longitudes:	                                                            
        alt, az = gal_to_topo(l, b)
        
        # make sure the alt, az is within the telescopes boundries
        if (leusch.ALT_MIN < alt < leusch.ALT_MAX) and \
           (leusch.AZ_MIN < az < leusch.AZ_MAX):      
            LT.point(alt, az)  
            spectrometer.read_spec(file_name + str(l) + ".fits", N, (l, b), 'ga')
				
	#take on and off noise data that we will use for temperature calibration at the end 
	capture_noise_calibration(file_name, N, l, b)
	LT.stow()
	#save the values of the longtitude were data is missing
	longitudes_values_missing(file_name)

def capture_missing_galactic_data(
    file_name,
    missing_data_file_name = "missing_longitudes_values",
    N=10
):
    '''	
    Capture long term missing data using the parameters longtidude range mentioned bellow
    '''	

	
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
