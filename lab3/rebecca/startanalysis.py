import ugradio
import numpy as np
import matplotlib.pyplot as plt
from ugradio import timing

#returns dictionary with arrays of collected data
def load_saves(filename):
    a = np.load(filename)
    return dict(zip(("{}".format(k) for k in a), (a[k] for k in a)))

def plotsun(time,data):
    '''Takes the raw output from the multimeter and gives back its plot and the Fourier Transform.

    Inputs:
    time = array of times from the multimeter
    data = array of voltages from the multimeter

    Outputs:
    A graph displaying the data and its transform
    freq = The frequencies from the transform, calculated using the time array
    power = The Fourier transform of the data array'''
    
    timefix,datafix= [],[]
    for i in range(0,1500):
        timefix.append(time[i])
        datafix.append(data[i])
    for j in range(1650,len(time)):
        timefix.append(time[j])
        datafix.append(data[j])

    timefix = np.array(timefix)
    datafix = np.array(datafix)
    
    trans = np.fft.fft(datafix)
    power = np.abs(trans)**2
    freq = np.fft.fftfreq(timefix.shape[-1])

    plt.subplot(1,2,1)
    plt.plot(timefix,datafix)
    plt.xlabel('Time (s)')
    plt.ylabel('Voltage (V)')

    plt.subplot(1,2,2)
    plt.semilogy(freq,np.fft.fftshift(power))
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Arbitrary Units (log(V$^2$))')

    plt.tight_layout()
    plt.show()

    return freq, power

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
    return A * np.cos(hour_angle) - B * np.sin(hour_angle)
