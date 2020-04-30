import ugradio 
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

def getline(file1,file2):
    '''Returns the callibrated spectrum for a pointing at a single LO frequency

    Inputs:
    file1: The name of the data file without noise as a string
    file2: The name of the data file with noise as a string

    Outputs:
    gainave * meanoff: The callibrated spectrum for one LO frequency.'''
    f_noise = fits.open(file1)
    f_quiet = fits.open(file2)
    
    on0,off0 = [],[]
    on1,off1 = [],[]

    on,off = [],[]
    for i in range(1,10):
        on0.append(f_noise[i].data['auto0_real'])
        off0.append(f_quiet[i].data['auto0_real'])

        on1.append(f_noise[i].data['auto1_real'])
        off1.append(f_quiet[i].data['auto1_real'])

        on.append(f_noise[i].data['auto0_real'])
        on.append(f_noise[i].data['auto1_real'])
        off.append(f_quiet[i].data['auto0_real'])
        off.append(f_quiet[i].data['auto1_real'])

    meanon0 = np.mean(on0,axis=0)
    meanoff0 = np.mean(off0,axis=0)

    meanon1 = np.mean(on1,axis=0)
    meanoff1 = np.mean(off1,axis=0)

    meanon = np.mean(on,axis=0)
    meanoff = np.mean(off,axis=0)

    gain0 = float(sum(meanoff0) * 80 / (sum(meanon0-meanoff0)))
    gain1 = float(sum(meanoff1) * 270 / (sum(meanon1-meanoff1)))
    gainave = .5*(gain0 + gain1)

    return gainave * meanoff

def filterline(line):
    '''Takes off a third degree polynomial from the baseline of the data

    Inputs:
    line: The spectrum needing to be filtered

    Outputs:
    line-baseline: The filtered spectrum'''
    x = np.array(range(len(line)))
    keepx = []
    keepline = []
    for i in range(0,5000):
        keepx.append(x[i])
        keepline.append(line[i])
    for k in range(6000,len(line)):
        keepx.append(x[i])
        keepline.append(line[i])
    base = np.polyfit(x,line,3)
    baseline = base[0]*x*x*x + base[1]*x*x + base[2]*x + base[3]
    return line-baseline

def getmeta(filename):
    '''Extracts the angle l, the right ascension, the declination, and the julian date of a data set

    Inputs:
    filename: The name of the data file as a string

    Outputs:
    l = Galactic longitude of data set
    ra = Right ascension of data set
    dec = Declination of data set
    jd = Julian date of data set'''
    f = fits.open(filename)
    dict(f[0].header)
    l =f[0].header[18]
    ra = f[0].header[20]
    dec = f[0].header[21]
    jd = f[0].header[22]
    return l, ra, dec, jd


def fanfan(label, start_angle, stop_angle):
    '''
    Will automatically plot every calibrated spectrum on the range.
    To see the next plot, close the existing one.
    This function would only really be helpful as a survey routine,
    i.e. giving the data a first glance.

    ex: fanfan('cycle_auto', 111, 250)
    you probably should not run this outside of a shell.
    '''
    for i in range(start_angle, stop_angle + 1):
        plt.plot(ffan(label, i)[4])
        plt.show()

def ffan(label, angle):
    '''
    Automatically expand a single reference into four file
    names for insertion into the function 'getfinal':
    ex ffan('cycle_auto', 250)
    '''
    n = '_noisy.fits'
    nn = '_quiet.fits'
    lo1 = '_634MHz'
    lo2 = '_635MHz'
    pre = label + '_' + str(float(angle)) + '_degrees'
    return getfinal(pre + lo1 + nn, pre + lo1 + n,
                    pre + lo2 + nn, pre + lo2 + n)

def getfinal(file1,file2,file3,file4):
    '''Returns a filtered spectrum that accounts for two LO frequencies

    Inputs:
    file1: The name of the data file without noise for the first LO frequency as a string
    file2: The name of the data file with noise for the first LO frequency as a string
    file3: The name of the data file without noise for the second LO frequency as a string
    file4: The name of the data file with noise for the second LO frequency as a string

    Outputs:
    l: Average galactic longitude between file1 and file3
    ra: Average right ascension between file1 and file3
    dec = Average declination between file1 and file3
    jd = Average julian date between file1 and file3
    final_line = The final spectrum after callibrations, filtering, and frequency switching'''
    line1 = getline(file1,file2)
    line2 = getline(file3,file4)
    filtered1 = filterline(line1)
    filtered2 = filterline(line2)
    final_line = filtered1-filtered2
    l1,ra1,dec1,jd1 = getmeta(file1)
    l2,ra2,dec2,jd2 = getmeta(file3)
    l = .5*(l1+l2)
    ra = .5*(ra1+ra2)
    dec = .5*(dec1+dec2)
    jd = .5*(jd1+jd2)
    return l,ra,dec,jd,final_line
