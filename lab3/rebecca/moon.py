import numpy as np
import ugradio
import matplotlib.pyplot as plt

#returns dictionary with arrays of collected data
def load_saves(filename):
    a = np.load(filename)
    return dict(zip(("{}".format(k) for k in a), (a[k] for k in a)))

def getmoon():
    dat = load_saves('moon_h2h_final.npz')
    data,stamp = dat['data'],dat['stamp']
    moov,moot = data[0],data[1]
    return moot,moov,stamp

def ha(hour,minute,second):
    return (hour * 3600 + minute * 60 + second) / 240

def deg(degrees, minutes, seconds, radians = False):
    return degrees + minutes/60 + seconds / 3600

#CURRENTLY ONLY CONSIDERING MOON

def hour_angle(day):
    '''Hard coded for Moon'''
    ra,dec = ugradio.coord.moonpos(day)
    hs = ugradio.timing.lst(day) - ra
    if hs < 0:
        hs += 360
    if hs > 360:
        hs += -360
    return hs, dec

def geodel(day):
    Bew = 20 #m
    lam = .025 #m
    hs,dec = hour_angle(day)
    return (Bew * np.cos(dec) / lam) * np.sin(hs)

def local_fringe(day):
    hs,dec = hour_angle(day)
    geod = geodel(day) / np.sin(hs)
    conversion = (2*np.pi)/(24*60*60)
    return conversion * geod * np.cos(hs)

def local_plot(stamp):
    time = []
    loc = []
    hs = []
    for i in range(len(stamp)):
        jd = ugradio.timing.julian_date(stamp[i][2])
        time.append(jd)
        h,dec = hour_angle(jd)
        hs.append(h)
        loc.append(local_fringe(jd))
    plt.plot(time,loc)
    plt.xlabel('Time (jd)')
    plt.ylabel('Local Fringe Frequency (Hz)')
    plt.show()
