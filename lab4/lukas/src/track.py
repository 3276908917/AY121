# for a backup that uses hard-coded values, see the latest edition in Rebecca's folder

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
        self.lat = ugradio.leo.lat
        self.lon = ugradio.leo.lon
        
        self.telescope = leusch.LeuschTelescope()
        self.noise = leusch.LeuschNoise()
        self.spec = leusch.Spectrometer()

        #I'm a little confused on how this one works since it has SynthDirect and then SynthClient
        self.lo = ugradio.agilent.SynthDirect()

    def find_point_safe(self, el, be):
        '''
        Calculate alt and az from galactic coordinates and
        check to make sure they are within bounds
        '''
        success = True
        now_jd = ugradio.timing.julian_date()
        alt, az = gal_to_topo(el, be, self.lat, self.lon, now_jd)
        if alt < leusch.ALT_MIN or alt > leusch.ALT_MAX or \
           az < leusch.AZ_MIN and az > leusch.AZ_MAX:
            print('Pointing out of bounds:')
            print('(l, b, alt, az) = ', el, be, alt, az)
            success = False
        return alt, az, success

    def single_measurement(self, el, be, full_prefix, N):
        alt_target, az_target, valid = self.find_point_safe(el, be)
        now = ugradio.timing.local_time()
        
        if valid:
            self.telescope.point(alt_target, az_target)
            alt_true, az_true = self.telescope.get_pointing()

            self.noise.on()

            self.lo.set_frequency(634, 'MHz')
            self.spec.read_spec(full_prefix + '_634MHz_noisy.fits', N, (el, be))
            
            self.lo.set_frequency(635, 'MHz')
            self.spec.read_spec(full_prefix + '_635MHz_noisy.fits', N, (el, be))

            self.noise.off()
            
            self.lo.set_frequency(634, 'MHz')
            self.spec.read_spec(full_prefix + '_634MHz_quiet.fits', N, (el, be))
            
            self.lo.set_frequency(635, 'MHz')
            self.spec.read_spec(full_prefix + '_635MHz_quiet.fits', N, (el, be))
            
        else:
            alt_true = az_true = None

        return np.array([el, be, alt_target, az_target, alt_true, az_true, now, valid])

    def scan_collect(self, list_targets, label, N=10):
        '''
        Collect @N spectra
        by observing each (galactic) coordinate pair in list_targets
            (currently handles only degrees)
        and save the data in two files, each with @label in the name
        the .fits file stores the spectra
        the .npz file stores the actual and desired pairs of topocentric coordinates,
            as well as the intended galactic latitude and current time.
        '''

        # It may be more dangerous to check the connection only once,
            # but heiles is slow and we want to reduce the number of unnecessary calculations
        self.spec.check_connected()
        meta_record = []
        
        for coordinate_pair in list_targets:
            el = coordinate_pair[0]
            be = coordinate_pair[1]
            meta_record.append(self.single_measurement(
                el, be, label + '_' + str(el) + '_degrees', N)
            )
            
        np.savez(label + '_stamp', stamp=meta_record)

        print('Ready. If you are done, remember to stow.')

        #self.telescope.stow()

    # possibly deprecated, transform_maps approach would be recommended
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
