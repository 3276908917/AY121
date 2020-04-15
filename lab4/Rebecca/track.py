import ugradio
import numpy as np

# I know there's a better way to do this, but, I don't know it at the moment
# So here's everything relevant from rotations.py
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

def gal_to_topo(el, be,
    lat=ugradio.nch.lat, lon=ugradio.timing.nch.lon,
    jd=ugradio.timing.julian_date(), radians=False
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
    lst = ugradio.timing.lst(jd, lon)
    hrd = np.dot(np.linalg.inv(M_eq_to_ha(lst)), ra_dec)
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

import ugradio.leusch as leusch

class Plane():
    def __init__(self):
        # temporary hard-coding
        self.lat = 37.91934
        self.lon = -122.15385
        
        self.telescope = leusch.LeuschTelescope()
        self.noise = leusch.LeuschNoise()
        self.spec = leusch.Spectrometer()

        #I'm a little confused on how this one works since it has SynthDirect and then SynthClient
        #self.lo = ugradio.agilent.SynthDirect()

    def find_point(self, el, be):
        '''
        Calculate alt and az from galactic coordinates and
        check to make sure they are within bounds
        '''
        alt, az = gal_to_topo(el, be, self.lat, self.lon)
        assert alt <= leusch.ALT_MIN and alt >= leusch.ALT_MAX and \
           az <= leusch.AZ_MIN and az >= leusch.AZ_MAX, \
            'Pointing out of bounds'
        return alt, az

    def pos_error(self, alt, az):
        '''
        Calculates the error between the wanted alt and az and
        the real alt and az of the telescope and returns those errors
        '''
        alt_real, az_real = self.telescope.get_pointing()
        return alt - alt_real, az - az_real

    def collect(self, el, be, full_prefix, N):
        list_alt_err, list_az_err = []
        
        for i in range(N):
            alt, az = self.find_point(el, be)
            self.telescope.point(alt, az)

            alt_err, az_err = self.pos_error(alt, az)
            list_alt_err.append(alt_err)
            list_az_err.append(az_err)

            assert self.spec.check_connected() == True, 'connection lost'
            ra, dec = gal_to_eq(el, be)
            self.spec.read_spec('plane_on_' + label + '_' +
                str(i) + '.fits', N, (ra, dec), 'eq')

        return list_alt_err, list_az_err

    def take_data(self, el, be, label, N=10):
        '''
        Collect @N spectra
        by observing galactic coordinates (@el, @be)
            (currently handles only degrees)
        and save the data in two files, each with @label in the name
        the .fits file stores the spectra
        the .npz file stores the two angle-errors for the pointings.
        '''
        self.lo.set_frequency(634, 'MHz')
        on_alt_err, on_az_err = self.collect(el, be, 'plane_off_' + label, N)

        self.lo.set_frequency(635, 'MHz')
        off_alt_err, off_az_err = self.collect(el, be, 'plane_on_' + label, N)

        np.savez('err_' + label,
                 on_alt_e=on_alt_err, on_az_e=on_az_err,
                 off_alt_e=off_alt_err, off_az_e=off_az_err)
        
        self.telescope.stow()

    def visibility_check(self, el, be, times, radians=False):
        '''
        Given a list of times in Julian format, return
        the times when the wanted galactic coordinates are
        within view of the telescope
        '''
        verified_times = []
            
        for t in times:
            alt, az = gal_to_topo(el, be, self.lat, self.lon, t, radians=False)            
            if alt >= leusch.ALT_MIN and alt <= leusch.ALT_MAX and \
               az >= leusch.AZ_MIN and az <= leusch.AZ_MAX:
                verified_times.append(times[i])

        return verified_times
