# off.npz was saved at 155850 on 2/25/20
# on.npz was saved at 152907 on 2/25/20

import numpy as np
import ugradio.timing

eq2g = np.array([[-.054876, -.873437, -.483835], [.494109, -.444830, .746982], [-.867666, -.198076, .455984]])

def eq_to_ha(LST):
    s = np.sin(LST)
    c = np.cos(LST)
    return np.array([[c, s, 0], [-s, c, 0], [0, 0, 1]])

#def topocentric(az, alt):
#    return np.array([np.cos(az) * np.cos(alt), np.sin(az) * np.cos(alt), np.sin(alt)])

def ha_to_topo(phi):
    s = np.sin(phi)
    c = np.cos(phi)
    return np.array([[-s, 0, c], [0, -1, 0], [c, 0, s]])

def rct_gal(l, b):
    return np.array([np.cos(b) * np.cos(l), np.cos(b) * np.sin(l), np.sin(b)])

def gal_to_topo(el, be, lat, radians=False):
    if not radians:
        l = np.radians(el)
        b = np.radians(be)
        phi = np.radians(lat)
    else:
        l = el
        b = be
        phi = np.radians(lat)
    rct = rct_gal(l, b)
    ra_dec = np.dot(np.linalg.inv(eq2g), rct)
    hrd = np.dot(np.linalg.inv(eq_to_ha(ugradio.timing.lst())), ra_dec)
    topo = np.dot(ha_to_topo(phi), hrd)
    return new_sphere(topo, radians)

def new_sphere(out_arr, radians=False):
    gp = np.arctan2(out_arr[1], out_arr[0])
    tp = np.arcsin(out_arr[2])
    if not radians:
        return np.degrees(gp), np.degrees(tp)   
    return gp, tp
