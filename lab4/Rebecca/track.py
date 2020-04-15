import ugradio
import numpy as np

#I know there's a better way to do this, but, I don't know it at the moment
#So here's everything relevant from rotations.py
M_eq_to_gal = np.array([
    [-.054876, -.873437, -.483835],
    [.494109, -.444830, .746982],
    [-.867666, -.198076, .455984]
])

def rectangle(a, b):
    '''
    Given a pair of angles (both angles must be in radians),
    return the corresponding 3x1 rectangular vector.
    '''
    return np.array([np.cos(b) * np.cos(a), np.cos(b) * np.sin(a), np.sin(b)])

def M_eq_to_ha(LST=ugradio.timing.lst()):
    '''
    Return the change-of-basis matrix between the equatorial and
    hour angle declination coordinate systems.
    The conversion depends on the @LST, Local Siderial Time
    '''
    s = np.sin(LST)
    c = np.cos(LST)
    return np.array([[c, s, 0], [s, -c, 0], [0, 0, 1]])

def M_ha_to_topo(phi=np.radians(ugradio.nch.lat)):
    '''
    Return the change-of-basis matrix between the hour angle declination
    and topocentric coordinate systems.
    The conversion depends on the user's current latitude @phi,
        which must be given in radians.
    '''
    s = np.sin(phi)
    c = np.cos(phi)
    return np.array([[-s, 0, c], [0, -1, 0], [c, 0, s]])

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

def gal_to_topo(el, be, lat=ugradio.nch.lat, radians=False):
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
        phi = lat
    rct = rectangle(l, b)
    ra_dec = np.dot(np.linalg.inv(M_eq_to_gal), rct)
    hrd = np.dot(np.linalg.inv(M_eq_to_ha(ugradio.timing.lst())), ra_dec)
    topo = np.dot(M_ha_to_topo(phi), hrd)
    return new_sphere(topo, radians)

def gal_to_eq(el, be, lat=ugradio.nch.lat, radians=False):
    if not radians:
        l = np.radians(el)
        b = np.radians(be)
        phi = np.radians(lat)
    else:
        l = el
        b = be
        phi = lat
    rct = rectangle(l, b)
    ra_dec = np.dot(np.linalg.inv(M_eq_to_gal), rct)
    return new_sphere(ra_dec, radians)

#in the end, the one we care about is gal_to_topo. End taking from rotiations.py

class plane():
    def __init__(self):
        self.telescope = ugradio.leusch.LeuschTelescope
        self.noise = ugradio.leusch.LeuschNoise()
        self.spec = ugradio.leusch.Spectrometer()

        #I'm a little confused on how this one works since it has SynthDirect and then SynthClient
        self.lo = ugradio.leusch.SynthDirect()
        self.lo.set_frequency(643, 'MHz')


    def pointing(el,be):
        '''Calculates alt and az from galactic coordinates and checks to make sure they are within bounds'''
        topo, rad = gal_to_topo(el,be)
        alt,az = topo[0],topo[1]
        if np.degrees(alt) <= ugradio.leusch.ALT_MIN or\
           np.degrees(alt) >= ugradio.leusch.ALT_MAX or\
           np.degrees(az) <= ugradio.leusch.AZ_MIN or\
           np.degrees(az) >= ugradio.leusch.AZ_MAX:
            print('Pointing out of bounds')
        elif:
            return alt,az

    def pointtime(el, be, lat=ugradio.nch.lat, radians=False,times):
        '''Given a list of times in julian_date, returns the times when the wanted galactic coordinates are within view of the telescope'''
        approved = []
        if not radians:
            l = np.radians(el)
            b = np.radians(be)
            phi = np.radians(lat)
        else:
            l = el
            b = be
            phi = lat
        rct = rectangle(l, b)
        for i in range(len(times)):
            ra_dec = np.dot(np.linalg.inv(M_eq_to_gal), rct)
            hrd = np.dot(np.linalg.inv(M_eq_to_ha(ugradio.timing.lst(times[i]))), ra_dec)
            topo = np.dot(M_ha_to_topo(phi), hrd)
            alt,az = new_sphere(topo, radians)
            if np.degrees(alt) >= ugradio.leusch.ALT_MIN and\
               np.degrees(alt) <=ugradio.leusch.ALT_MAX and\
               np.degrees(az) >= ugradio.leusch.AZ_MIN and\
               np.degrees(az) <= ugradio.leusch.AZ_MAX:
                approved.append(times[i])

    def position(alt,az):
        '''Positions the telecope given an alt and az in radians'''
        alt,az = np.degrees(alt),np.degrees(az)
        self.telescope.point(alt,az)

    def pos_error(alt,az):
        '''Calculates the error between the wanted alt and az and the real alt and az of the telescope and returns those errors'''
        alt,az = np.degrees(alt),np.degrees(az)
        alt_real,az_real = self.telescope.get_pointing()
        return alt-alt_real,az-az_real

    def take_date(el,be,label,N=10):
        alt_err, az_err = []
        count = 0
        for i in range(N):
            while count < N:
                alt, az = pointing(el, be)
                position(alt, az)
                alt_err,az_err = pos_error(alt, az)
                alterr.append(alt_err)
                azerr.append(az_err)
                if self.spec.check_connected() == True:
                    radec, rad = gal_to_eq(el, be)
                    self.spec.read_spec('planeon' + label + str(count) + '.fits', N, radec, 'J2000')
                count += 1

        self.lo.set_frequency(644, 'MHz')
        count = 0
        for j in range(N):
            while count < N:
                alt, az = pointing(el, be)
                position(alt, az)
                alt_err, az_err = pos_error(alt,az)
                alterr.append(alt_err)
                azerr.append(az_err)
                if self.spec.check_connected() == True:
                    radec, rad = gal_to_eq(el, be)
                    self.spec.read_spec('planeoff' + label + str(count) + '.fits', N, radec, 'J2000')
                count += 1
        self.telescope.stow()
