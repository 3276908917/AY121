import ugradio 
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

#the goal of this script is to examine the contributions and limitations of
#the two different polarizations on the leuschner telescope
#My conlcusion: one polarization is much cleaner, but the other signal is stronger
#This whill change depenidng on the pointing, but one of the two will always
#Have the stronger signal. I do recommend using both as the final result is still pretty clean

#files are already hardcoded in, so make sure that the following datafiles
#are in the same directory as the script:
#thirty_plane_on.fits
#thirty_plane_off.fits
#attempt_all_140.0_degrees_634MHz_noisy.fits
#attempt_all_140.0_degrees_634MHz_quiet.fits
#attempt_all_140.0_degrees_635MHz_noisy.fits
#attempt_all_140.0_degrees_635MHz_quiet.fits

def polargraphs():
    '''Compares the polarization responses in plots

    Outputs:
    Top Left Graph: The averaged spectra from auto0_real
    Top Right Graph: The averaged spectra from auto1_real
    Bottom Left Graph: The averages spectra when auto0_real and auto1_Real are added linearly
    Bottom Right Graph: An overlapyed graoh of the previous three spectra'''
    fon = fits.open('thirty_plane_on.fits')
    foff = fits.open('thirty_plane_off.fits')
    mon,moff = [],[]
    mon0,moff0 = [],[]
    mon1,moff1 = [],[]
    for i in range(1,10):
        mon.append(fon[i].data['auto0_real'])
        mon.append(fon[i].data['auto1_real'])
        
        mon0.append(fon[i].data['auto0_real'])
        mon1.append(fon[i].data['auto1_real'])
        
        moff.append(foff[i].data['auto0_real'])
        moff.append(foff[i].data['auto1_real'])
        
        moff0.append(foff[i].data['auto0_real'])
        moff1.append(foff[i].data['auto1_real'])
    meanon = np.mean(mon,axis=0)
    meanon0 = np.mean(mon0,axis=0)
    meanon1 = np.mean(mon1,axis=0)
    
    meanoff = np.mean(moff,axis=0)
    meanoff0 = np.mean(moff0,axis=0)
    meanoff1 = np.mean(moff1,axis=0)
    
    sline = meanon/meanoff
    sline0 = meanon0/meanoff0
    sline1 = meanon1/meanoff1

    plt.subplot(2,2,1)
    plt.plot(sline0,label='auto0_real')
    plt.ylabel('Arbitrary Units Proportional to Power')
    plt.legend()
    
    plt.subplot(2,2,2)
    plt.plot(sline1,label='auto1_real')
    plt.legend()
    
    plt.subplot(2,2,3)
    plt.plot(sline,label='auto0_real + auto1_real')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Arbitrary Units Proportional to Power')
    plt.legend()

    plt.subplot(2,2,4)
    plt.plot(sline0,label='auto0_real')
    plt.plot(sline1,label='auto1_real')
    plt.plot(sline,label='auto0_real + auto1_real')
    plt.xlabel('Frequency (Hz)')
    plt.legend()

    plt.show()

def lo1():
    '''Compares the gains of the two polarizations for the first LO value

    Outputs:
    A graph that overlays the averages spectra of both polarizations with gain, and then the final averaged spectra with both polarizations added linearly
    np.array(gain0*meanoff0): callibrated spectra with auto0_real
    np.array(gain1*meanoff1): callibrated spectra with auto1_real
    np.array(gainave*meanoff): callibrated spectra with both polarizations'''
    
    f_noise = fits.open('attempt_all_140.0_degrees_634MHz_noisy.fits')
    f_quiet = fits.open('attempt_all_140.0_degrees_634MHz_quiet.fits')
    
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

    plt.plot(gain0*meanoff0,label='auto0 only')
    plt.plot(gain1*meanoff1,label='auto1 only')
    plt.plot(gainave*meanoff,label='auto0+auto1')
    plt.legend()
    plt.show()
    return np.array(gain0*meanoff0),np.array(gain1*meanoff1),np.array(gainave*meanoff)

def lo1_filter():
    '''This applies a third degree polynomial to the data for the first lo and subtracts it off in order ot eliminate the largest features of the systematics

    Outputs:
    Graph that displays the filtering with auto0_real,auto1_real, and then both polarizations added linearly
    line0-baseline0: The filtered spectra for auto0_real
    line1-baseline1: The filtered spectra for auto1_real
    line-baseline: The filtered spectra for both polarizations'''
    line0,line1,line = lo1()
    x = np.array(range(len(line)))
    base0 = np.polyfit(x,line0,3)
    base1 = np.polyfit(x,line1,3)
    base = np.polyfit(x,line,3)
    baseline0 = base0[0]*x*x*x + base0[1]*x*x + base0[2]*x + base0[3]
    baseline1 = base1[0]*x*x*x + base1[1]*x*x + base1[2]*x + base0[3]
    baseline = base[0]*x*x*x + base[1]*x*x + base[2]*x + base[3]
    #plt.plot(line,label="Unfiltered, auto0_real+auto1_real")
    plt.plot(line0-baseline0,label="auto0_real")
    plt.plot(line1-baseline1,label="auto1_real")
    plt.plot(line-baseline,label="auto0_real+auto1_real")
    plt.xlabel('Index')
    plt.ylabel('Proportional to Temperature (K)')
    plt.legend()
    plt.show()
    return line0-baseline0,line1-baseline1,line-baseline

def lo2():
    '''Compares the gains of the two polarizations for the second LO value

    Outputs:
    A graph that overlays the averages spectra of both polarizations with gain, and then the final averaged spectra with both polarizations added linearly
    np.array(gain0*meanoff0): callibrated spectra with auto0_real
    np.array(gain1*meanoff1): callibrated spectra with auto1_real
    np.array(gainave*meanoff): callibrated spectra with both polarizations'''
    f_noise = fits.open('attempt_all_140.0_degrees_635MHz_noisy.fits')
    f_quiet = fits.open('attempt_all_140.0_degrees_635MHz_quiet.fits')
    
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

    plt.plot(gain0*meanoff0,label='auto0 only')
    plt.plot(gain1*meanoff1,label='auto1 only')
    plt.plot(gainave*meanoff,label='auto0+auto1')
    plt.legend()
    plt.show()
    return np.array(gain0*meanoff0),np.array(gain1*meanoff1),np.array(gainave*meanoff)

def lo2_filter():
    '''This applies a third degree polynomial to the data for the second lo and subtracts it off in order ot eliminate the largest features of the systematics

    Outputs:
    Graph that displays the filtering with auto0_real,auto1_real, and then both polarizations added linearly
    line0-baseline0: The filtered spectra for auto0_real
    line1-baseline1: The filtered spectra for auto1_real
    line-baseline: The filtered spectra for both polarizations'''
    line0,line1,line = lo2()
    x = np.array(range(len(line)))
    base0 = np.polyfit(x,line0,3)
    base1 = np.polyfit(x,line1,3)
    base = np.polyfit(x,line,3)
    baseline0 = base0[0]*x*x*x + base0[1]*x*x + base0[2]*x + base0[3]
    baseline1 = base1[0]*x*x*x + base1[1]*x*x + base1[2]*x + base0[3]
    baseline = base[0]*x*x*x + base[1]*x*x + base[2]*x + base[3]
    plt.plot(line0-baseline0,label="auto0_real")
    plt.plot(line1-baseline1,label="auto1_real")
    plt.plot(line-baseline,label="auto0_real+auto1_real")
    plt.xlabel('Index')
    plt.ylabel('Proportional to Temperature (K)')
    plt.legend()
    plt.show()
    return line0-baseline0,line1-baseline1,line-baseline

def outcome():
    '''Showd final outcomes from filtering. Warning: You will get all the graphs from previous functions unless you comment them out'''
    fit00,fit10,fit0 = lo1_filter()
    fit01,fit11,fit1 = lo2_filter()
    plt.plot(fit00-fit01,label="auto0_real")
    plt.plot(fit10-fit11,label="auto1_real")
    plt.plot(fit0-fit1,label="auto0_real+aut1_real")
    plt.legend()
    plt.show()
