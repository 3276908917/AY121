import numpy as np
import ugradio
import matplotlib.pyplot as plt

#returns dictionary with arrays of collected data
def load_saves(filename):
    a = np.load(filename)
    return dict(zip(("{}".format(k) for k in a), (a[k] for k in a)))

def getsun():
    '''Loads one hour sun trial.

    Inputs:
    NA

    Outputs:
    sut = array of times
    suv = array of voltages'''
    dat = load_saves('sun_63.73_minutes.npz')
    data = dat['data']
    suv,sut = data[0],data[1]
    return sut,suv

#from rotations.py
def ha(hour,minute,second):
    return (hour * 3600 + minute * 60 + second) / 240

#from rotations.py
def deg(degrees, minutes, seconds, radians = False):
    return degrees + minutes/60 + seconds / 3600

#CURRENTLY ONLY CONSIDERING SUN

def hour_angle(time):
    '''Derives the hour angle of the sun on March 12, 2020 for inputted time in unix time

    Inputs:
    time = a single value of time in unix time

    Outputs:
    hs = the hour angle of the sun at that time on March 12, 2020'''
    convert =  (15/1) 
    midday = 2458921.29236
    hs = convert * ((ugradio.timing.julian_date(time)-midday) * 24)
    return hs

def declination():
    '''Derives the declination of the sun based on how many days have past since Jan 1 of the current year. Harded for March 12th, 2020.

    Outputs:
    dec = The declination of the sun on the give day'''
    day = 71
    coeff = (360/365.24) * (np.pi/180)
    theta1 = -23.44 * (np.pi/180)
    theta2 = coeff * (day+10) + (2*.0167)*np.sin(coeff*(day-2))
    dec = np.arcsin((np.sin(theta1)*np.cos(theta2)))*(180/np.pi)
    return dec

def local_fringe(time):
    '''Calculates the local fringe value using Taylor expanded format and initial parameters provided by lab manual.

    Inputs:
    time = a single value of time in unix time

    Outputs:
    The value of the local fringe at the given time'''
    Bew = 20 #m
    lam = .025 #m
    hs = hour_angle(time) * (np.pi/180)
    dec = (declination()) * (np.pi/180)
    convert = (2*np.pi) / (24*60*60)
    return (convert * Bew *np.cos(dec) * np.cos(hs))/ (lam)

def local_plot(time):
    '''Plots local fringe values over a period of time.

    Inputs:
    time = array of times

    Outputs:
    A graph showing the local fringe frequency values over the given timespan'''
    times = []
    loc = []
    for i in range(len(time)):
        times.append((time[i]-1583996400)/(60*60))
        loc.append(local_fringe(time[i]))
    plt.plot(times,loc)
    plt.show()

#good enough for the moment. It has the right shape and the right spread, it's just off by about .01 Hz from the data, which seems pretty good.

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

def fftfilter(t,v,fil):
    '''Filters out the lowest 30 frequencies in the transform

    Inputs:
    t = array of times
    v = array of voltages
    fil = the number of frequencies to be filtered out

    Outputs:
    freq = array of frequencies
    trans = transform of voltages with lowest 30 frequencies filtered out'''
    freq,power,trans = getfft(t,v)
    for i in range(0,fil):
        trans[i] = 0
        trans[-i] = 0
    return freq, trans

def filterplot(t,v):
    '''Plots the filtered signal as well as the power spectrum of the filtered signal.

    Inputs:
    t = array of times
    v = array of voltages

    Outputs:
    Graph showing filtered signal and filtered power spectrum'''
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

def dualplot(t,v):
    '''Plots data and its power spectrum.
    Inputs:
    t = array of times
    v = array of voltages

    Outputs:
    Graph with signal and its power spectrum'''
    freq,power,trans = getfft(t,v)
    plt.subplot(2,1,1)
    plt.plot(t,v)
    plt.xlabel('Time (s)')
    plt.ylabel('Voltage (V)')

    plt.subplot(2,1,2)
    plt.plot(np.fft.fftshift(freq),np.fft.fftshift(power))
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Arbitrary Unites (log(V$^2$))')

    plt.tight_layout()
    plt.show()

#trying for that optimized fit
def geodel(time):
    '''Calculates the value of the modified geometric delay from the source to the interferometer. Hard Coded for March 12, 2020

    Inputs:
    time = a single value of time in unix time

    Outputs:
    The value of the modified geometric delay'''
    bew = 20 #m
    lam = .025 #m
    dec = declination()*(np.pi/180)
    hs = hour_angle(time)*(np.pi/180)
    return (bew * np.cos(dec) / lam) * np.sin(hs)

def fringe_data(t,v):
    timing = 600
    times = []
    fringes = []
    breaks = int(len(t)/timing)
    for i in range(breaks):
        time = t[i*timing: (i+1)*timing]
        times.append(t[i*timing])
        volts = v[i*timing:(i+1)*timing]
        freq,trans = fftfilter(time,volts,5)
        power = np.abs(trans)**2
        fringe = np.abs(freq[np.argmax(power)])
        fringes.append(fringe)
    return times,fringes



def brute(t,v,A,B):
    '''Calculates the sum of the risiduals of the data and a guess for Qew and Qns

    Inputs:
    x = array of times
    y = array of voltages
    A = guess for A
    B = guess for B

    Outputs:
    risid = array of the risuduals
    sum_tot = sum of the total of the squares of the rididual'''
    risid = []
    times,fringes = fringe_data(t,v)
    for i in range(len(times)):
        risid.append((fringes[i] - (A*np.cos(hour_angle(times[i])*(np.pi/180)) + B*np.sin(hour_angle(times[i])*(np.pi/180))))**2)
    sum_sq_tot = sum(risid)
    return risid, sum_sq_tot

def risidplot(t,v,A,B):
    '''Plots the filtered data, the fitted guess, and then the risidual

    Inputs:
    x = array of times
    y = array of voltages
    A = guess for A value
    B = guess for B value

    Outputs:
    Graph with data, fit, and risidual'''
    risid = []
    times,fringes = fringe_data(t,v)
    date = []
    for i in range(len(times)):
        date.append(12 + hour_angle(times[i])/15)
        risid.append((fringes[i] - (A*np.cos(hour_angle(times[i])*(np.pi/180)) + B*np.sin(hour_angle(times[i])*(np.pi/180))))**2)
    plt.subplot(2,1,1)
    plt.plot(date,fringes,label="Data")
    plt.plot(date,local_fringe(times),label="Initial Values")
    plt.plot(date,A*np.cos(hour_angle(times)*(np.pi/180)) + B*np.sin(hour_angle(times)*(np.pi/180)),label="Optimized Values")
    plt.ylabel('Frequency (Hz)')
    plt.legend()
    plt.subplot(2,1,2)
    plt.plot(date,risid)
    plt.xlabel('Time of Day (hr)')
    plt.ylabel('Risidual (Hz)')
    
    plt.show()


def optim(t,v):
    '''Attempts to optimize a fit to the data

    Inputs:
    x = array of times
    y = array of voltages

    Outputs:
    Graph of Qew vs the sum of the square of the risidual'''
    #Bew = 20 #m
    #lam = .025 #m
    #dec = (declination()-.2) * (np.pi/180)
    #convert = (2*np.pi) / (24*60*60)
    #A = convert * Bew * np.cos(dec) / lam
    A = 0.03016748195476793 #these are the correct values
    B = 0.0001900000000000002
    dA = .0001
    dB = .00001
    As = []
    Bs = []
    sums = []
    #A = A - 10*dA
    Bnew = B - 10*dB
    for i in range(20):
        #As.append(A)
        Bs.append(Bnew)
        ___,trial = brute(t,v,A,Bnew)
        sums.append(trial)
        #A = A + dA
        Bnew = Bnew + dB
    print(Bnew[np.argmin(sums)])
    print(Bnew)
    plt.plot(Bnew,sums)
    plt.show()
    return A,B

def getB(A,B):
    lam = .025
    convert = (2*np.pi)/(24*60*60)
    dec = declination() * (np.pi/180)
    Bew = A * lam / (np.cos(dec) * convert)
    
    L = 37.873199 * (np.pi/180)
    Bns = B * lam / (np.cos(dec) * np.sin(L) * convert)
    return Bew,Bns
