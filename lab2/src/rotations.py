# off.npz was saved at 155850 on 2/25/20
# on.npz was saved at 152907 on 2/25/20

import numpy as np
import ugradio.timing

# The change-of-basis matrix between equatorial and galactic coordinate systems
eq_to_gal = np.array([[-.054876, -.873437, -.483835], [.494109, -.444830, .746982], [-.867666, -.198076, .455984]])

def eq_to_ha(LST):
    '''
    Return the change-of-basis matrix between the equatorial and
    hour angle declination coordinate systems.
    The conversion depends on the @LST, Local Siderial Time
    '''
    s = np.sin(LST)
    c = np.cos(LST)
    return np.array([[c, s, 0], [s, -c, 0], [0, 0, 1]])

def ha_to_topo(phi):
    '''
    Return the change-of-basis matrix between the hour angle declination
    and topocentric coordinate systems.
    The conversion depends on the user's current latitude @phi,
        which must be given in radians.
    '''
    s = np.sin(phi)
    c = np.cos(phi)
    return np.array([[-s, 0, c], [0, -1, 0], [c, 0, s]])

def rct_gal(l, b):
    '''
    Given a pair of angles @l, @b (both angles must be in radians),
    return the corresponding 3x1 rectangular vector.
    '''
    return np.array([np.cos(b) * np.cos(l), np.cos(b) * np.sin(l), np.sin(b)])

def gal_to_topo(el, be, lat, radians=False):
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
    else:
        l = el
        b = be
        phi = np.radians(lat)
    rct = rct_gal(l, b)
    ra_dec = np.dot(np.linalg.inv(eq_to_gal), rct)
    # The program at least works fine up until here
    hrd = np.dot(np.linalg.inv(eq_to_ha(ugradio.timing.lst())), ra_dec)
    topo = np.dot(ha_to_topo(phi), hrd)
    # new_sphere has also been demonstrated to work
    return new_sphere(topo, radians)

def new_sphere(out_arr, radians=False):
    '''
    Given a 3x1 vector,
    return the corresponding pair of angles
    @radians determines whether the angles are given in radians.
    '''
    gp = np.arctan2(out_arr[1], out_arr[0])
    tp = np.arcsin(out_arr[2])
    if not radians:
        return np.degrees(gp), np.degrees(tp)   
    return gp, tp
