# galactic plane is at 0 degrees (galactic longitude)

import numpy as np
import ugradio.leo as loc
import ugradio.leusch as leusch
import matplotlib.pyplot as plt

plane = 0
LST = np.linspace(0, 360, num = 360)

def bound_altaz(bound):
    return [gal_to_topo_lst(bound, plane, lst,
        loc.lat, loc.lon, radians = False) \
            for lst in LST]
        
# modified form from rotations.py
def gal_to_topo_lst(el, be, lst, lat, lon, radians=False
):
    '''
    @radians determines the format of BOTH input and output!
    Given a pair of angles @el and @be (in galactic coordinates),
    return a pair of angles relating the associated
    azimuth and altitude.
    '''
    if not radians:
        l = np.radians(el)
        b = np.radians(be)
        phi = np.radians(lat)
        theta = lon
    else:
        l = el
        b = be
        phi = lat
        theta = np.degrees(lon)
    rct = rectangle(l, b)
    ra_dec = np.dot(np.linalg.inv(M_eq_to_gal), rct)
    hrd = np.dot(np.linalg.inv(M_eq_to_ha(lst)), ra_dec)
    topo = np.dot(M_ha_to_topo(phi), hrd)
    return new_sphere(topo, radians)

def bound_plot(ell):
    fixed_ell = np.array(bound_altaz(ell))
    alts = fixed_ell[:, 0]
    azs = fixed_ell[:, 1]
    
    fig, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=True)

    plt.xlabel('LST [degrees]', fontsize=12)
    fig.text(0, 0.5, r'Angle [degrees]',
             va='center', rotation='vertical', fontsize=12)

    ax1.tick_params(axis="x", labelsize=12)
    ax2.tick_params(axis="x", labelsize=12)
    ax1.tick_params(axis="y", labelsize=12)
    ax2.tick_params(axis="y", labelsize=12)

    ax1.plot(LST, alts, label='altitude')
    ax1.plot(LST, [leusch.ALT_MIN] * len(LST), label='minimum allowed')
    ax1.plot(LST, [leusch.ALT_MAX] * len(LST), label='maximum allowed')
    ax1.legend(loc='upper right')
    
    ax2.plot(LST, azs, label='azimuth')  
    ax2.plot(LST, [leusch.AZ_MIN] * len(LST), label='minimum allowed')
    ax2.plot(LST, [leusch.AZ_MAX] * len(LST), label='maximum allowed')
    ax2.legend(loc='upper right')
    
ELL = np.linspace(-10, 250, 260)

def time_altaz(lst):
    return [gal_to_topo_lst(ell, plane, lst,
        loc.lat, loc.lon, radians = False) \
            for ell in ELL]

def gal_topo_plot(lst):
    fixed_lst = np.array(time_altaz(lst))
    alts = fixed_lst[:, 0]
    azs = fixed_lst[:, 1]
    
    fig, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=True)

    plt.xlabel('LST [degrees]', fontsize=12)
    fig.text(0, 0.5, r'Angle [degrees]',
             va='center', rotation='vertical', fontsize=12)

    ax1.tick_params(axis="x", labelsize=12)
    ax2.tick_params(axis="x", labelsize=12)
    ax1.tick_params(axis="y", labelsize=12)
    ax2.tick_params(axis="y", labelsize=12)

    ax1.plot(ELL, alts, label='altitude')
    ax1.plot(ELL, [leusch.ALT_MIN] * len(ELL), label='minimum allowed')
    ax1.plot(ELL, [leusch.ALT_MAX] * len(ELL), label='maximum allowed')
    ax1.legend(loc='upper right')
    
    ax2.plot(ELL, azs, label='azimuth')  
    ax2.plot(ELL, [leusch.AZ_MIN] * len(ELL), label='minimum allowed')
    ax2.plot(ELL, [leusch.AZ_MAX] * len(ELL), label='maximum allowed')
    ax2.legend(loc='upper right')
