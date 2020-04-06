#!/usr/bin/env python
# coding: utf-8

# In[12]:


import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import dft
from scipy import fft
import scipy.stats as stats
import scipy
import math
import ugradio
from ugradio import timing
from itertools import zip_longest
import matplotlib.gridspec as gridspec

def mean(signal):
    return np.mean(signal, axis=0)

def median(signal):
    return np.median(signal, axis=0)

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

def inverse_fourier(signal):
    return np.fft.ifft(signal)

# returns dictionary with arrays of collected data
def load_saves(filename):
    a = np.load(filename)
    return dict(zip(("{}".format(k) for k in a), (a[k] for k in a)))

def Histogram(Signal, bins=15):

    mu = mean(Signal)
    sigma= mean(Signal)

    # print(mu, sigma)

    # data = total_rms
    #plt.figure(figsize=(10,10))
    plt.subplot(221)
    # Plot the histogram.
    plt.hist(Signal, bins, density=True, alpha=0.6, color='b')
    # Fit a normal distribution to the data:
    mu, std = norm.fit(Signal)
    # Plot the PDF.
    xmin, xmax = plt.xlim()
    x = np.linspace(xmin, xmax, 100)
    p = norm.pdf(x, mu, std)
    plt.plot(x, p, 'k', linewidth=2, label='Bins = 15')
    # title = "mean = %.2f,   std = %.2f" % (mu, std)
    # plt.title(title)

    plt.ylabel('Nomalized Distribution')
    plt.legend()
    plt.show()

    return x


# In[13]:


data_moon = load_saves('moon_part2_final.npz')

data_moon_final = load_saves('moon_h2h_final.npz')

data_sun_800 = load_saves('sundata_892.73_minutes.npz')

data_CygA = load_saves('CygA_21.45_minutes.npz')

data_moon_75 = load_saves('moon_h2h_437.12_minutes.npz')

data_CygA_final = load_saves('CygA_final.npz')

data_sun = load_saves('sun_63.73_minutes.npz')

data_final = load_saves('sh3_final.npz')

data_sun_746 = load_saves('sun_h2h_746.1_minutes.npz')


# In[15]:



def ha(hour,minute,second):
    #credit to Rebecca gore 
    return (hour * 3600 + minute * 60 + second) / 240

def deg(degrees, minutes, seconds, radians = False):
    #credit to Rebecca gorn 
    return degrees + minutes/60 + seconds / 3600

def hour_angle(time):
    
     #credit to Rebecca gore
        
    ''' converts unix time to hour angle '''
    
    hourangle = 15 * ((ugradio.timing.julian_date(time)- 2458921.29236) * 24)
    return hourangle

def declination(day=71):
    
    #credit to Rebecca gore 
    ''' Get declination in degrees from julian date'''
    
    degree_to_radian = (np.pi/180)
    radian_to_degree = (180/np.pi)
    
    coeff = (360/365.24) * degree_to_radian
    
    theta1 = -23.44 * degree_to_radian
    
    theta2 = coeff * (day + 10) + (2*.0167)*np.sin(coeff*(day-2))
    
    dec = np.arcsin((np.sin(theta1)*np.cos(theta2)))
    
    return dec * radian_to_degree

def Calculated_local_fringe_plot(time):
    #credit goes to Rebecca gore 
    '''Plots local fringe values over a period of time.
    Inputs:
    time = array of times
    Outputs:
    A graph showing the local fringe frequency values over the given timespan'''
    times = []
    loc = []
    for i in range(len(time)):
        times.append((time[i]- 1583996400)/(60*60))
        loc.append(local_fringe_frequency(time[i])-.015)
    plt.plot(times,loc)
    plt.show()


# In[159]:



def Plot_real_imaginary(signal, time):

    N = len(signal) # Sample count

    t = time # Time vector


    # FFT + bins + normalization
    bins = np.fft.fftfreq(N, st)    
    fft  = [i / (N/2) for i in fourier(signal)]

    return bins, np.real(fft) ,np.imag(fft)

def local_fringe_frequency(time):
       
    #get the local fringe frequency in cycles per seconds
    
    degree_to_radian = (np.pi/180)
    radian_to_degree = (180/np.pi)

    Bew = 2000 #m
    Bns = 0

    lam = .025 #m

    hs = hour_angle(time) * degree_to_radian

    dec = (declination(71)-.2) * degree_to_radian

    seconds_converter = (2*np.pi) / (24*60*60)

    Local_fringe = (( seconds_converter * Bew / (lam)) * np.cos(dec)) * np.cos(hs)


    return Local_fringe


def fourier_frequency(signal, time):
    
    #calculate power spectrum
    
    sig = np.fft.fft(signal)
    power = np.abs(sig)**2
    freq = np.fft.fftfreq(time.shape[-1])
    
    return freq, sig


def power_frequency(signal, time):
    
    #calculate power spectrum
    
    sig = np.fft.fft(signal)
    power = np.abs(sig)**2
    freq = np.fft.fftfreq(time.shape[-1])
    
    return freq, power

def Filter(signal, time, filter_frequency = 8):
    
    #filter the signal from the frequencies -filter_frequency to +filter_frequency
    freq, sig =  fourier_frequency(signal, time)
    
    for i in range(0,filter_frequency):
        sig[i], sig[-i] = 0
              
    return freq, sig

def Plot_filtered_signal(signal, time, filter_frequency=8 ):
   
    freq, sig = Filter(signal, time, filter_frequency)
    
    Filtered_signal = inverse_fourier(sig)
            
    fig, (ax1, ax2) = plt.subplots(1, 2)

    #Plot signal after filtering 
    ax1.plot(time , Filtered_signal)
    plt.xlabel('Time (s)')
    plt.ylabel('Voltage (V)')
    
    #plot fft transform of filtered signal
    power = np.abs(sig)**2
    ax2.plot(np.fft.fftshift(freq), np.fft.fftshift(power))
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Voltage squared (V$^2$)')

    plt.show()

def plot_data_and_fft(time,data):

    time = np.array(time)
    data = np.array(data)

    sig_fft = fourier(data)
    power = np.abs(sig_fft)**2
    
    freq = np.fft.fftfreq(time.shape[-1])
    
    fig, (ax1, ax2) = plt.subplots(1, 2)

    ax1.plot(time, data)
    plt.xlabel('Time (s)')
    plt.ylabel('Voltage (V)')
    
    ax2.plot(np.fft.fftshift(freq),np.fft.fftshift(power))
    plt.yscale('log')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Amplitude (log(V$^2$))')
    


    return freq, power

#path delay calculated from baseline, declination, hour angle and wavelenght 
def Path_delay(time):
    
    degree_to_radian = (np.pi/180)
    
    bew = 2000 #cm
    lam = 2.5 #cm
    
    dec = declination()*degree_to_radian
    
    return (bew * np.cos(dec) / lam) 


# In[160]:


print(declination())
# print (hour_angle)


# In[161]:


#Amplitude [0] and time [1]
Data_moon_75 = data_moon_75['data']

Data_moon = data_moon['data']

Data_moon_final = data_moon_final['data']

Data_sun_800 = data_sun_800['data']

Data_final = data_final['data']

Data_CygA = data_CygA['data']

Data_CygA_final = data_CygA_final['data']

Data_sun_746 = data_sun_800['data']

Data_sun = data_sun['data']


#Stamp of the data (Altitude, Azimuth and unix time)

stamp_moon_75 = data_moon_75['stamp']

stamp_moon_final = data_moon_final['stamp']

stamp_sun = data_sun_800['stamp']

stamp_sun_800 = data_sun_800['stamp']

stamp_CygA = data_CygA['stamp']

stamp_CygA_final = data_CygA_final['stamp']

stamp_CygA_final_stamp = data_CygA_final['stamp']

stamp_sun_746 = data_sun_800['stamp']


# In[162]:



# plot_data(Data_final[1], Data_final[0])


# In[163]:



plt.plot(np.cos(hour_angle(stamp_sun[2])*(np.pi/180)))
# plt.plot(np.sin(2*np.pi*hour_angle(stamp_sun[2])))
# plt.plot()
# plt.plot(hour_angle(stamp_sun[2]))


# In[164]:


# from the lab jupiter nootbook
def plt_fringe(bl):
    plt.figure(figsize=(8,8))
    plt.clf()
    ax2 = plt.subplot(111,projection='polar')
    fringe_pattern = np.exp(-2j * np.pi * bl * np.sin(theta) / c * fq)
    plt.plot(theta, 0.1 * fringe_pattern.real + 1)
    ax2.set_theta_offset(np.pi/2)
    ax2.set_yticks([])
    ax2.set_thetamin(-90)
    ax2.set_thetamax(90)
    plt.title('Fringe versus Angle')
    plt.show()

    
def plt_beam(point):
    plt.figure(figsize=(9,9))
    plt.clf()
    ax2 = plt.subplot(111,projection='polar')
    beam_pattern = np.exp(-(theta+point)**2 / (2 * 0.1**2)) # Gaussian beam pattern
    fringe_pattern = np.exp(-2j * np.pi * bl * np.sin(theta) / c * fq)
    plt.polar(theta, 0.1 * fringe_pattern.real + 1)
    plt.polar(theta, 0.1 * beam_pattern * fringe_pattern.real + 1)
    plt.polar(theta, beam_pattern)
    ax2.set_theta_offset(np.pi/2)
    ax2.set_yticks([])
    ax2.set_thetamin(-90)
    ax2.set_thetamax(90)
    plt.title('Response versus Angle')
    plt.show()
    
from ipywidgets import interact

SIZE = 2000

theta = np.linspace(-np.pi/2, np.pi/2, SIZE)

Bew = 2000 #cm
Bns = 0 #cm
lam = 2.5 #cm

Second_converter = 2 * np.pi / (60 * 60 * 24)

bl = 2000. # baseline, cm east-west
fq = 10.7e9 # spectral frequency, Hz
c = 3e10 # speed of light, cm / s




i = interact(plt_fringe, bl=(500, 3000))


# In[165]:


i = interact(plt_beam, point=(-np.pi/2, np.pi/2))


# In[166]:


get_ipython().run_line_magic('matplotlib', 'notebook')

def Qew(value, Range, number_of_samples):
    
    '''Creates guesses of Qew values around the first guessed value were the baseline(ew)= 20m and baseline(ns)=0'''
    x= value-(Range/2)
    return np.linspace(x, x+Range, number_of_samples)

def Qns(value, Range, number_of_samples):
    ''' Creates guesses of Qns values around the first a given value '''
    
    x= value-(Range/2)
    return np.linspace(x, x+Range, number_of_samples)



def Baseline_ew(Qew, declination = declination(), wavelength=0.025):
    return ((Qew * wavelength)/ np.cos(declination))


def Baseline_ns(Qns, declination=declination(), wavelength=0.025, terrestrial_latitude=37.873199):
    return (Qns * wavelength / np.cos(declination) / np.sin(np.radians(terrestrial_latitude)))
#     return (Qns * wavelength / np.cos(declination) / np.sin(np.radians(terrestrial_latitude)))* ((60 * 60 * 24) / (2 * np.pi))


def baseline_value(Bew, Bns):
    
    return np.sqrt(Bew**2 + Bns**2)

#path delay calculated from baseline, declination, hour angle and wavelenght 
def Path_delay(hs, Qew, Qns):
    
    return Qew * np.sin(hs) + Qns * np.cos(hs)


def least_square_A_B(signal, time, Qew, Qns):
    '''    
    F(hs) = A*X + B*Y
    Y=sin(2*pi*v*t)
    X=cos(2*pi*v*t)
    
    '''
    hour_angles = np.array(np.radians(hour_angle(time)))

    path_delay= Path_delay(hour_angles, Qew, Qns)



    X, Y = np.cos(path_delay), np.sin(path_delay)
    
    sig = np.sum(signal * hour_angles)
    
   
    A = np.vstack([X, Y , np.ones(len(X))]).T
    
#     print ('A', len(A))
#     print ('X', len(X))
#     print ('Y', len(Y))
#     print ('sig', len(signal))

    lsq = np.linalg.lstsq(A, signal, rcond=None)
#     lsq = np.linalg.lstsq(A, sig, rcond=None)
    print(lsq)

    A, B = lsq[0][0], lsq[0][1]

    residual = lsq[1]

    minimum_value_cordinate = lsq[2]

    return A, B

def sum_of_squares(signal, time, Qew, Qns):
    """
    Calculating the residual squared from the guessed A and B, its the sum of ( signal - (Acos(2*pi*path_delay) + Bsin(2*pi*path_delay)))^2
    """
    hs = np.array(np.radians(hour_angle(time)))
    
    A, B = least_square_A_B(signal, time, Qew, Qns)

    path_delay = Path_delay(hs, Qew, Qns)
    
    residuals = signal - (A*np.cos(path_delay) + B*np.sin(path_delay))
    
    return np.sum(residuals**2)


def baseline_1Dimension(signal, time, declination=declination()):
    
#   Q_ew =  2000 / 2.5 (cos(-3.326465607576737) * (2*pi/60*60*24)

    Q_e_w = 798.6520965
    Qew_guessed_values = np.linspace(770,830, 3824) 
    
    Qns_guessed_values, Bns = 0, 0
    
    residual_squared = [sum_of_squares(signal, time, Q, Qew_guessed_values) for Q in Qew_guessed_values]
    
    min_Qew, min_s_squared = minimum_1D(Qew_guessed_values, residual_squared)

    Bew = Baseline_ew(min_Qew, declination)
    
    baseline = baseline_value(Bew, Bns)

    print("Qew value: " , str(min_Qew))
    print("S_squared value: " , str(min_s_squared))
    print("Baseline: " , str(baseline))
    

    plt.plot(Qew_guessed_values, residual_squared, color='b', label= ' Q_ew value vs residual squared' )
    plt.scatter(min_Qew, min_s_squared, color='r', label= 'minimum value')
    
    plt.xlabel(" Q_ew Values")
    plt.ylabel("Residual square Values S$^2$")
    plt.legend()
  
    plt.show()
    
    return residual_squared, min_Qew, min_s_squared


def baseline_2Dimensions(signal , time, declination=declination()):

#     Qew_guessed_values = Qew(0.05807961942, .1, 2000) 
#     Qns_guessed_values = Qns(0.1, 1, 2000)

    Q_e_w = 798.6520965



    Qew_guessed_values = np.linspace(770,830, 3824)
    Qns_guessed_values = np.linspace(600, 900, 3824) 
    residual_squared = []
    
    for Q__ns in Qns_guessed_values:
        Qns_dimension = []
        for Q__ew in Qew_guessed_values:
            Qns_dimension.append(sum_of_squares(signal, time, Q__ew, Q__ns))
            
        residual_squared.append(Qns_dimension)
        
    residual_squared = np.array(residual_squared)
    
    min_Qew, min_Qns, min_s_squared = minimum_2D(Qew_guessed_values, Qns_guessed_values, residual_squared)
    Bew, Bns = Baseline_ew(min_Qew, declination), Baseline_ns(min_Qns, declination)
    baseline = baseline_value(Bew, Bns)
    
    print("Qew value: "  , str(min_Qew), "Qns value: " , str(min_Qns))
    print("S_squared value: " , str(min_s_squared))
    print("Baseline: " , str(baseline))
    print ('B_ew and B_ns', Bew, Bns)
    
    # 2 d meshcolor the plot
    fig, ax = plt.subplots()
    cmap = ax.pcolormesh(Qew_guessed_values, Qns_guessed_values, residual_squared)
    fig.colorbar(cmap)
    plt.show(fig)
    
    return Bew, Bns, baseline





# In[167]:


print (Data_sun[0])
baseline_1Dimension(Data_sun[1], Data_sun[0], declination=declination())


# In[143]:


plt.plot(Data_sun[1], Data_sun[0])
plt.show()


# In[ ]:


baseline_2Dimensions(Data_sun[0], Data_sun[1], declination=declination())


# In[ ]:





# In[ ]:




