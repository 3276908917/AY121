import ugradio
import numpy as np
import matplotlib.pyplot as plt
from ugradio import timing

def fourier_skeleton(x, y, xBounds=None, yBounds=None, logv=False,
    xLabel='Frequency (Hz)', yLabel = r'Voltage (mV)', rude_filter = False):
    fig = plt.figure(figsize=(6,3))
    plt.subplots_adjust(left=.15, bottom=.15, right=.95, top=.9)

    ax = fig.add_subplot(111)

    ax.tick_params(axis="x", labelsize=12)
    ax.tick_params(axis="y", labelsize=12)

    if rude_filter:
        y[0:1] = 0

    ax.plot(x, np.real(np.fft.fftshift(y * 1000)), label='real')
    ax.plot(x, np.imag(np.fft.fftshift(y * 1000)), label='imaginary')
    
    plt.xlabel(xLabel, fontsize=12)
    plt.ylabel(yLabel, fontsize=12)
    if xBounds is not None:
        plt.xlim(xBounds)
    if yBounds is not None:
        plt.ylim(yBounds)
    if logv:
        plt.yscale('log')
    plt.show()

    # Now we want to find the locations of the peaks
    print(np.abs(x[np.argmax(np.fft.fftshift(y))]))

# hard-code city
def collect_fringes(volts):
    fringes = []
    # We evaluate by ten minute segments (10 * 60 = 600 seconds per)
    x = fr(1, 600)
    P = lambda r : np.abs(np.fft.fft(r)) ** 2
    for i in range(600, len(volts), 600):
        y = volts[i - 600:i]
        y = P(y)
        y[0:1] = 0 # naive filtering approach
        y = np.fft.fftshift(y)
        fringes.append(np.abs(x[np.argmax(y)]))
    #print(fringes)
    return fringes

def pp3(x, y, xBounds=None, yBounds=None, logv=False, xLabel='Frequency (MHz)',
    yLabel = r'Magnitude-squared Voltage (V$^2$)'):
    '''
    Plot @y versus @x where @y is shifted to center the zero frequency
        (more at numpy.fft.fftshift)
    The plot is automatically labeled such that the
        x-axis corresponds to frequency in megahertz
        y-axis corrosponds to magnitude-squared voltage in square volts
    One can include tuples @xBounds and @yBounds to zoom the graph appropriately
        where each tuple is in the form (lower_bound, upper_bound).
    One can ahead-of-time specify the y-axis to be logarithmic
        with @logv=True. However, it is not clear whether this functionality is useful,
        because one can simply press L (while focused on a matplotlib plot)
        to achieve the same effect.
    '''
    fig = plt.figure(figsize=(6,3))
    plt.subplots_adjust(left=.15, bottom=.15, right=.95, top=.9)

    ax = fig.add_subplot(111)

    ax.tick_params(axis="x", labelsize=12)
    ax.tick_params(axis="y", labelsize=12)

    ax.plot(x, np.fft.fftshift(y))
    plt.xlabel(xLabel, fontsize=12)
    plt.ylabel(yLabel, fontsize=12)
    if xBounds is not None:
        plt.xlim(xBounds)
    if yBounds is not None:
        plt.ylim(yBounds)
    if logv:
        plt.yscale('log')
    plt.show()

'''

Notes to self:
    sunburst seems to stop at about 1.5853214e9
        and resumes at about        1.58535819e9
    indices 11k to 42k seem safe

    since we took data at 1 Hz, the recommended 10 minutes
        will correspond to 600 indices.

'''

''' uNtr HIs puNkt -ai wrk aser> b RbekO '''

#returns dictionary with arrays of collected data
def load_saves(filename):
    a = np.load(filename)
    return dict(zip(("{}".format(k) for k in a), (a[k] for k in a)))

def freq_range(v_s, N, W=1):
    """
    Return a N-length array
    Where frequencies range between plus or minus W * v_s / 2
    """
    lobe = round(N / 2)
    interval = W * v_s / N
    return np.array([i * interval for i in range(-lobe, lobe)])

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
    print(timefix.shape[-1])
    freq = np.fft.fftfreq(timefix.shape[-1])
    print(freq)

    plt.subplot(1,2,1)
    plt.plot(timefix,datafix)
    plt.xlabel('Time (s)')
    plt.ylabel('Voltage (V)')

    plt.subplot(1,2,2)
    x = freq_range(1, len(time) - 150)
    
    plt.semilogy(np.fft.fftshift(freq), np.fft.fftshift(power))
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
