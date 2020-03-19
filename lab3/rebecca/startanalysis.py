import ugradio
import numpy as np
import matplotlib.pyplot as plt
from ugradio import timing

#returns dictionary with arrays of collected data
def load_saves(filename):
    a = np.load(filename)
    return dict(zip(("{}".format(k) for k in a), (a[k] for k in a)))

def getfft(t,v):
    '''Calculates the fourier transform of a voltage data set.
    Inputs:
    t = array-like of times
    v = array-like of voltages

    Outputs:
    freq = array-like of the corresponding frequencies
    power = arra-like of the fourier transform'''
    trans = np.fft.fft(v)
    power = np.abs(trans)**2
    freq = np.fft.fftfreq(t.shape[-1])
    return freq, power

def plotsun(time,data):
    '''Takes the raw output from the multimeter and gives back its plot and the Fourier Transform. Edit for initial sun trial hardcoded.

    Inputs:
    time = array of times from the multimeter
    data = array of voltages from the multimeter

    Outputs:
    A graph displaying the data and its transform
    freq = The frequencies from the transform, calculated using the time array
    power = The Fourier transform of the data array'''
    
    timefix,datafix= [],[] #hardcoded for 1 hour data
    for i in range(0,1500):
        timefix.append(time[i])
        datafix.append(data[i])
    for j in range(1501,1649):
        timefix.append(0)
        datafix.append(0)
    for k in range(1650,len(time)):
        timefix.append(time[k])
        datafix.append(data[k])

    timefix = np.array(timefix)
    datafix = np.array(datafix)
    
    freq, power = gettfft(timefix,datafix)

    plt.subplot(2,1,1)
    plt.plot(timefix,datafix)
    plt.xlabel('Time (s)')
    plt.ylabel('Voltage (V)')

    plt.subplot(2,1,2)
    plt.plot(np.fft.fftshift(freq),np.fft.fftshift(power))
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Arbitrary Units (log(V$^2$))')

    plt.tight_layout()
    plt.show()

    return freq, power

def getplot(t,v):
    time = []
    for i in range(len(t)):
        jd = ugradio.timing.julian_date(t[i])
        time.append(ugradio.timing.lst(jd))
    plt.plot(time,v)
    #plt.xlabel('Time (s)')
    #plt.ylabel('Voltage (V)')

def dualplot(t,v):
    freq,power = getfft(t,v)
    plt.subplot(1,2,1)
    getplot(t,v)

    plt.subplot(1,2,2)
    plt.semilogy(np.fft.fftshift(freq),np.fft.fftshift(power))
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Arbitrary Unites (log(V$^2$))')

    plt.tight_layout()
    plt.show()

def fourplot(t1,t2,v1,v2):
    '''For plotting both initial moon trials. Edits hardcoded'''
    t1,v1 = t1[:3000],v1[:3000]
    t2,v2 = t2[1000:],v2[1000:]
    freq1,power1 = getfft(t1,v1)
    freq2,power2 = getfft(t2,v2)
    
    plt.subplot(2,2,1)
    getplot(t1,v1)
    plt.xlabel('Time (s)')
    plt.ylabel('Voltage (V)')

    plt.subplot(2,2,2)
    getplot(t2,v2)
    plt.xlabel('Time (s)')

    plt.subplot(2,2,3)
    plt.semilogy(np.fft.fftshift(freq1),np.fft.fftshift(power1))
    plt.xlabel('First Set; Frequency (Hz)')
    plt.ylabel('Arbitrary Unites (log(V$^2$))')

    plt.subplot(2,2,4)
    plt.semilogy(np.fft.fftshift(freq2),np.fft.fftshift(power2))
    plt.xlabel('Second Set; Frequency (Hz)')

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
