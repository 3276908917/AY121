import ugradio
import numpy as np
import matplotlib.pyplot as plt
from ugradio import timing

#dc = ugradio.interf_delay.DelayClient()
#dc.delay_ns(0.)

#returns dictionary with arrays of collected data
def load_saves(filename):
    a = np.load(filename)
    return dict(zip(("{}".format(k) for k in a), (a[k] for k in a)))

def getsun():
    sun = load_saves('suntrial_772.16_minutes.npz')
    sundat = sun['data']
    sunv,sunt = sundat[0],sundat[1]
    sunv,sunt = sunv[1750:41000],sunt[1750:41000]
    return sunt, sunv

def getfft(t,v):
    '''Calculates the fourier transform of a voltage data set.
    Inputs:
    t = array of times
    v = array of voltages

    Outputs:
    freq = array of the corresponding frequencies
    power = array of the fourier transform'''
    trans = np.fft.fft(v)
    power = np.abs(trans)**2
    freq = np.fft.fftfreq(t.shape[-1])
    return freq, power, trans

def fftfilter(t,v):
    '''Filters out lowest thirty frequencies of the transform, positive and negative.
    
    Inputs:
    t = array of times
    v = array of voltages
    
    Outputs:
    freq = array of frequencies
    trans = array of filtered transform values
    '''
    freq,power,trans = getfft(t,v)
    for i in range(0,30):
        trans[i] = 0
        trans[-i] = 0
    return freq, trans

def filterplot(t,v):
    '''Plots the filtered signal and its transform
    
    Inputs:
    t = array of times
    v = array of voltages
    
    Outputs:
    Graph of both the filtered signal and its transform    
    '''
    freq,trans = fftfilter(t,v)
    filtered = np.fft.ifft(trans)
    power = np.abs(trans)**2

    plt.subplot(2,1,1)
    plt.plot(t,filtered)
    plt.xlabel('Time (s)')
    plt.ylabel('Voltage (V)')

    plt.subplot(2,1,2)
    plt.plot(np.fft.fftshift(freq),np.fft.fftshift(power))
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Arbitrary Unites (log(V$^2$))')

    plt.show()

def getplot(t,v):
    #time = []
    #for i in range(len(t)):
    #    jd = ugradio.timing.julian_date(t[i])
    #    time.append(ugradio.timing.lst(jd))
    plt.plot(t,v)
    plt.xlabel('Time (s)')
    plt.ylabel('Voltage (V)')

def dualplot(t,v):
    freq,power = getfft(t,v)
    plt.subplot(1,2,1)
    getplot(t,v)

    plt.subplot(1,2,2)
    plt.plot(np.fft.fftshift(freq),np.fft.fftshift(power))
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Arbitrary Unites (log(V$^2$))')

    plt.tight_layout()
    plt.show()

def declination(day):
    '''Derives the declination of the sun based on how many days have past since Jan 1 of the current year

    Inputs:
    day = The number of days since Jan 1 of the current year

    Outputs:
    dec = The declination of the sun on the give day
    dec_real = dec_real is the actual declination of the sun
    dec-dec_real = gives the difference between the calculation and real value'''
    jd = ugradio.timing.julian_date()
    ra,dec_real = ugradio.coord.sunpos(jd)
    coeff = (360/365.24) * (np.pi/180)
    theta1 = -23.44 * (np.pi/180)
    theta2 = coeff * (day+10) + (2*.0167)*np.sin(coeff*(day-2))
    dec = np.arcsin((np.sin(theta1)*np.cos(theta2)))*(180/np.pi)
    return dec,dec_real, dec-dec_real

def chisq(y,yexp,sigma):
    return np.sum(np.abs(y-yexp)**2/sig**2)
                  

def brute():
    Qew = (Bew/lam) * np.cos(delta)
    Qns = (Bns/lam) * np.sin(L) * np.cos(delta)
    

def thefringe():
    geodel = 0 #place holder
    theta = 2 * np.pi * geodel
    fringe_amp = A * np.cos(theta) + B * np.sin(theta)

def getlocal(A,B):
    Bew = 20 #m
    Bns = 0 #m
    lam = .025 #
    ra,dec = ugradio.sunpos(ugradio.timing.julian_date())
    hour_angle = ugradio.timing.lst() - ra
    return A * np.cos(hour_angle)
