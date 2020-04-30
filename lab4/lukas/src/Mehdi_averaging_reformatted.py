import ugradio
import scipy.signal 
from astropy.io import fits
import numpy as np
from scipy import fft
import scipy.stats as stats
import scipy
import matplotlib.pyplot as plt

'''
def rms(signal):
    return np.sqrt(mean(signal)**2)
'''

# This has to be redundant.
def power_spectrum(signal):
    fourier_transform= np.fft.fft(signal)
    return np.abs(fourier_transform)**2

def data_header(file_name):
    """ Get header from file name """
    leushner_header = fits.open(str(file_name))
    return leushner_header[0].header
		
def plot_raw_leusch(file_name):	
    """
    Graph animmated plots of both auto0_real and
    auto1_real over all the samples N
    """
    f = fits.open(str(file_name))
    N = f[0].header['NSPEC']

    for p in np.arange(1, N):
        plt.ion()

        plt.plot(f[p].data['auto0_real'], label=' auto0_real')
        plt.plot(f[p].data['auto1_real'], label= 'auto1_real')

        plt.xlabel('samples')
        plt.ylabel('Amplitude ($V^2$)')
        # plt.title(' Power spectrum ')
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
            # print ('leushner_data[i].data[name]', leushner_data[i].data[name])
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
    data = fits.open(str(file_name))
			
    header = data_header(file_name)

    N = header['NSPEC']
    
    #average of data, either first or second 
    average_spectra = Average_spectrum(data, N, polarization_number)
    
    print ('Average_spectra',average_spectra)
    
    #frequency responce of our AVERAGE data spectrum
    freq = Frequency_spectrum(header, average_spectra)
    
    plt.plot(freq, average_spectra, label=file_name)
    plt.xlabel('Frequency')
    plt.ylabel('$V^2$')
    plt.legend()
    # plt.show()

				
def Frequency_spectrum(header, spectrum_data):
    """
    Changes samples axis of the spectrum data
    into frequency axis	
    """		
    vsamp = header["SAMPRATE"]
            
    t = 1/vsamp
    
    freq = (np.fft.fftshift(np.fft.fftfreq(len(spectrum_data), t ))) / 1e6 + 1420  # where 1420 is the total of the LO frequencies in MHz that shifted  down our HI signal, I added this to re-shift the received signal toward frequency 1420.4 mhz

    return freq


def base_fit_average(average_spectra, frequency):
    """
    Filter our average data, then use 1d interpolation
    on the data within boundaries 
    """	
    filtered_spectra = scipy.signal.medfilt(average_spectra, kernel_size=9) 
    max_index = np.argmax(filtered_spectra)
    
    bandwidth = 200
    mask = np.ones(len(filtered_spectra), dtype=bool)
    
    for i in range(len(mask)):
            if  (max_index - bandwidth) < i < (max_index + bandwidth) :
                    mask[i] = False
                    
    baseline_interpolation = scipy.interpolate.interp1d(frequency[mask], filtered_spectra[mask])
    
    return baseline_interpolation(frequency)

def doppler_velocity_correction(header, frequency, velocity_correction='off'):
    """
    Calculate the correction needed to shift the doppler velocity to the
    LSR (Local Standard of Rest) frame, which follows the mean motion of
    material in the Milky Way in the neighborhood of the Sun.
	
    We use the header from the collected data each time to get the ra, dec and jd;
    then we use them to get the doppler correction, then we get the doppler velocity.
    """
    # physical constants
    HI_freq = 1420.40575  # frequency of our HI signal in MHz
    c = 3e5  #speed of light in km/s 
		
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
        
        Velocity =  np.array(doppler_velocity) - (Doppler_correction.value/1e3) #convertion Doppler correction from m/s to km/s
	
    #turn off Doppler correction
    elif velocity_correction=='off':
        Velocity =  np.array(doppler_velocity)
             
    return  Velocity

def Gain(signal_noise_on, signal_noise_off, number_spectra=10, Noise_temp=301, Sky_temp=1):	
    """
    Calculate the Gain (G) like in Lab 2, uses the signal with and without noise, temperature of the sky and temprature from the noise in kelvin
    """
    noise_on_spectra = Average_spectrum(signal_noise_on, number_spectra)
    noise_off_spectra = Average_spectrum(signal_noise_off, number_spectra)
    # Noise_temp - Sky_temp = 300 #kelvin
    gain = ((Noise_temp - Sky_temp) / np.sum(noise_on_spectra - noise_off_spectra)) * np.sum(noise_off_spectra)
    return gain

def final_calibration(file_name, Polarization_number='first', Gain=1, velocity_correction='off' ):	
    """
    Calculate the gain first Gain(noise_on, noise_off, Noise_temp, Sky_temp)
    then use the final value of the gain for the calibration.
	
    Calibrate the leushner data, change it from frequency Mhz vs amplitude V^2 to temperature spectra (kelvin) with respect to doppler velocity m/s (LSR frame)
    """	
    gain = 300 # will calculate it later once	
    data = fits.open(str(file_name))
    header = data_header(file_name)
    N = header['NSPEC']	
	
    # Average of data, either first or second polarization. 
    average_spectra = Average_spectrum(data, N, Polarization_number)
    print ('Average_spectra', average_spectra)
    #frequency responce of our AVERAGE data spectrum
    freq = Frequency_spectrum(header, average_spectra)
    Temperature = (average_spectra - base_fit_average(average_spectra, freq)) * Gain
    print ('Temperature', Temperature)
    doppler_velocity = doppler_velocity_correction(header, freq, velocity_correction)
    return doppler_velocity, Temperature 
	
fon = fits.open('attempt_all_150.0_degrees_635MHz_noisy.fits')
N = fon[0].header['NSPEC']
	
foff= fits.open('attempt_all_150.0_degrees_635MHz_quiet.fits')
N = foff[0].header['NSPEC']

print (foff[0].header)
gain = Gain(fon, foff, N, 301, 1)
#calculate the gain:
print ('Gain', Gain(fon, foff, N, 301, 1))

def plot_temperature_doppler(data, velocity_correction= 'off'):
    data_x = data[0]
    data_y = data[1]			
	
    plt.plot(data_x, data_y, label=' Doppler correction is '+ velocity_correction)
	
    plt.xlabel('Doppler velocity (km/s)')
    plt.ylabel('Temperature (kelvin)')
    plt.legend()

def plot_all_attempts(lo_freq1, lo_freq2, file_name='attempt_all', velocity_correction='off', Gain= 2000 ):
    degree_spacing = 2/2
    lontitude_min= -10
    longtitude_max = 250
    Longtitude_Spacing = (longtitude_max - (lontitude_min)) / degree_spacing
	
    longitude_range = np.linspace(lontitude_min, longtitude_max, int(Longtitude_Spacing)+1) 
    print ('longitude_range',len(longitude_range))
		
    Missing_longitudes_values = np.array([])			
	
    for l in longitude_range:
        plt.ion()	
        try:		
            Lo1_noisy_name = file_name + '_' + str(float(int(l))) +'_degrees_'+ str(lo_freq1) +'MHz_noisy.fits'	
            Lo1_quiet_name = file_name + '_' + str(float(int(l))) +'_degrees_'+ str(lo_freq1) +'MHz_quiet.fits'
			
            Lo2_noisy_name = file_name + '_' + str(float(int(l))) +'_degrees_'+ str(lo_freq2) +'MHz_noisy.fits'	
            Lo2_quiet_name = file_name + '_' + str(float(int(l))) +'_degrees_'+ str(lo_freq2) +'MHz_quiet.fits'		
                                                        
            print(Lo1_noisy_name)	
            print(Lo1_quiet_name)

            Lo1_on = fits.open(Lo1_noisy_name)
            Lo1_off = fits.open(Lo1_quiet_name)

            Lo2_on = fits.open(Lo2_noisy_name)
            Lo2_off = fits.open(Lo2_quiet_name))

            gain_lo1 = 300 / np.sum(Average_spectrum(Lo1_on, N=10, polarization_number='first') \
                - Average_spectrum(Lo1_off, N=10, polarization_number='first')) * \
                np.sum(Average_spectrum(Lo1_off, N=10, polarization_number='first'))

            gain_lo2 = 300 / np.sum(Average_spectrum(Lo1_on, N=10, polarization_number='first') \
                - Average_spectrum(Lo1_off, N=10, polarization_number='first')) * \
                np.sum(Average_spectrum(Lo1_off, N=10, polarization_number='first'))
                        
            N = Lo1_on[0].header['NSPEC']
            print (N)

            data_lo1_on = final_calibration(Lo1_noisy_name, Polarization_number='first', Gain= gain_lo1, velocity_correction='off' )
            data_lo1_off = final_calibration(Lo1_quiet_name, Polarization_number='first', Gain= gain_lo1, velocity_correction='off' )

            data_lo2_on = final_calibration(Lo2_noisy_name, Polarization_number='first', Gain= gain_lo2, velocity_correction='off' )
            data2_lo2_off = final_calibration(Lo2_quiet_name, Polarization_number='first', Gain= gain_lo2, velocity_correction='off' )
                                
            #Plot average plots of all data in the file within the longtitude we collected
            #			plot_average(Lo1_noisy_name, polarization_number= 'first')

            #Plot Temperatude vs Doppler velocity of all data in the file within the longtitude we collected
            plot_temperature_doppler(data_lo1_on, velocity_correction= 'off')
            plot_temperature_doppler(data_lo2_on, velocity_correction= 'off')
                        
            plt.title('Longtitude'+ str(l))

            plt.pause(.1)
            plt.draw()
            plt.clf()
					
            except:
                print('Data is not taken for this longtitude')
                Missing_longitudes_values= np.append(Missing_longitudes_values, float(l))
	
    plt.ioff()	
    np.save("missing_longitudes_values", Missing_longitudes_values)
				
plot_all_attempts(634, 635, file_name = 'attempt_all')	

#plot_average('attempt_all_145.0_degrees_634MHz_noisy.fits', polarization_number= 'first' )
#plt.show()	

#data = final_calibration('attempt_all_160.0_degrees_635MHz_quiet.fits', velocity_correction='on', Gain= gain)
#plt.plot(data[0], data[1], label=' Doppler corrected' )
#plt.xlabel('Doppler velocity (km/s)')
#plt.ylabel('Temperature (kelvin)')
#plt.legend()
#plt.show()
#
