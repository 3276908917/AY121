# galactic plane is at 0 degrees (galactic longitude)

import numpy as np
import ugradio.leo as loc
import ugradio.leusch as leusch
import matplotlib.pyplot as plt
import rotations

plane = 0
list_LST = np.linspace(0, 360, num = 360)

def bound_altaz(bound_ell):
    '''
    Return a list of topocentric coordinate pairs
    each of which is a conversion of the galactic coordinates
    (ell = bound_ell, be = plane)
    at a different local sidereal time.
    The list covers one full sidereal day.
    '''
    return [gal_to_topo_lst(bound_ell, plane, lst,
        loc.lat, radians = False) \
            for lst in list_LST]
        
# modified form from rotations.py
def gal_to_topo_lst(ell, be, lst, lat, radians=False):
    '''
    Given a pair of angles @el and @be (in galactic coordinates),
    return a pair of angles relating the associated
    azimuth and altitude.

    @radians determines the format of BOTH input and output!
    '''
    if not radians:
        l = np.radians(ell)
        b = np.radians(be)
        phi = np.radians(lat)
        sidereal_angle = np.radians(lst)
    else:
        l = ell
        b = be
        phi = lat
        sidereal_angle = lst
    rct = rotations.rectangle(l, b)
    ra_dec = np.dot(np.linalg.inv(rotations.M_eq_to_gal), rct)
    # LST must be in radians
    hrd = np.dot(np.linalg.inv(rotations.M_eq_to_ha(sidereal_angle)), ra_dec)
    topo = np.dot(rotations.M_ha_to_topo(phi), hrd)
    raw_angles = rotations.new_sphere(topo, radians)
    alt = raw_angles[1]
    az = raw_angles[0]
    if az < 0:
        az += 360
    return az, alt

def bound_plot(bound_ell):
    '''
    Plot the topocentric coordinate conversions for all
    local sidereal times in a sidereal day.
    We hold the galactic longitude 
    For convenience, we also plot the maximum and minimum
    permitted values for the Leuschner dish pointing routine.
    '''
    top_over_time = np.array(bound_altaz(bound_ell))
    alts = top_over_time[:, 1]
    azs = top_over_time[:, 0]
    
    fig, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=True)

    plt.xlabel('LST [degrees]', fontsize=12)
    fig.text(0, 0.5, r'Angle [degrees]',
             va='center', rotation='vertical', fontsize=12)

    ax1.tick_params(axis="x", labelsize=12)
    ax2.tick_params(axis="x", labelsize=12)
    ax1.tick_params(axis="y", labelsize=12)
    ax2.tick_params(axis="y", labelsize=12)

    y_span = len(lst_LST)

    ax1.plot(list_LST, alts, label='altitude')
    ax1.plot(list_LST, [leusch.ALT_MIN] * y_span, label='minimum allowed')
    ax1.plot(list_LST, [leusch.ALT_MAX] * y_span, label='maximum allowed')
    ax1.legend(loc='upper right')
    
    ax2.plot(list_LST, azs, label='azimuth')  
    ax2.plot(list_LST, [leusch.AZ_MIN] * y_span, label='minimum allowed')
    ax2.plot(list_LST, [leusch.AZ_MAX] * y_span, label='maximum allowed')
    ax2.legend(loc='upper right')
    
list_ELL = np.linspace(-10, 250, 260)

def time_altaz(lst):
    '''
    Return a list of topocentric coordinate pairs
    each of which is a conversion of galactic coordinates.
    We hold the local sidereal time to be constant (given by @lst)
    as well as the galactic latitude (be = plane).
    We exclusively vary galactic longitude ell over the range
    of interest for this particular lab.
    '''
    return [gal_to_topo_lst(ell, plane, lst,
        loc.lat, radians = False) \
            for ell in list_ELL]

def gal_topo_plot(lst, anchor1='upper right', anchor2='upper right'):
    '''
    Plot the topocentric coordinate conversions for all
    galactic longitudes of interest
    at the given local sidereal time @lst.
    For convenience, we also plot the maximum and minimum
    permitted values for the Leuschner dish pointing routine.
    '''
    top_over_ell = np.array(time_altaz(lst))
    alts = top_over_ell[:, 1]
    azs = top_over_ell[:, 0]
    
    fig, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=False)

    plt.xlabel('Galactic Longitude [degrees]', fontsize=12)
    fig.text(0, 0.5, 'Topocentric Angle [degrees]',
             va='center', rotation='vertical', fontsize=12)

    ax1.tick_params(axis="x", labelsize=12)
    ax2.tick_params(axis="x", labelsize=12)
    ax1.tick_params(axis="y", labelsize=12)
    ax2.tick_params(axis="y", labelsize=12)

    ax1.plot(list_ELL, alts, label='altitude')
    ax2.plot(list_ELL, azs, label='azimuth')  

    y_span = len(list_ELL)
    
    ax1.plot(list_ELL, [leusch.ALT_MIN] * y_span, label='minimum allowed')
    ax1.plot(list_ELL, [leusch.ALT_MAX] * y_span, label='maximum allowed')
    ax2.plot(list_ELL, [leusch.AZ_MIN] * y_span, label='minimum allowed')
    ax2.plot(list_ELL, [leusch.AZ_MAX] * y_span, label='maximum allowed')

    ax1.legend(loc=anchor1)
    ax2.legend(loc=anchor2)

# deg_lst = np.degrees(ugradio.timing.lst(ugradio.timing.julian_date(ugradio.timing.unix_time() + 3600 * 0))) 
# gal_topo_plot(deg_lst); plt.show()
