import numpy as np
import ugradio
import ugradio.leusch as leusch

# make sure that you copy rotations.py into the same directory
import rotations

class Plane():
    def __init__(self, latitude=ugradio.leo.lat, longitude=ugradio.leo.lon):
        self.lat = latitude
        self.lon = longitude
        
        self.telescope = leusch.LeuschTelescope()
        self.noise = leusch.LeuschNoise()
        self.spec = leusch.Spectrometer()

        self.lo = ugradio.agilent.SynthDirect()

    def find_point_safe(self, el, be, quiet=False):
        '''
        Calculate alt and az from galactic coordinates and
        check to make sure they are within bounds.
        Instead of erroring-out,
        this version simply records and prints the command failure.
        '''
        success = True
        now_jd = ugradio.timing.julian_date()
        alt, az = rotations.gal_to_topo(
            el, be, now_jd, self.lat, self.lon, radians=False
        )
        if alt < leusch.ALT_MIN or alt > leusch.ALT_MAX or \
           az < leusch.AZ_MIN or az > leusch.AZ_MAX:
            if not quiet:
                print('Pointing out of bounds:')
                print('(l, b, alt, az) = ', el, be, alt, az)
            success = False
        return alt, az, success

    def single_measurement(self, el, be, full_prefix, N):
        '''
        Collect @N spectra
        by observing the galactic coordinates (@el, @be)
        which are first converted into topocentric coordinates.
        We generate four .fits files, for which we turn noise
        on and off, and switch the LO between 634 and 635 MHz.
        '''
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

        return np.array([alt_target, az_target, alt_true, az_true, valid])

    def scan_collect(self, list_targets, label, N=10):
        '''
        Collect @N spectra
        by observing each (galactic) coordinate pair in list_targets
            (currently handles only degrees)
        and save the data in files named according to @label
        the single .npz file stores the actual and desired pairs of topocentric coordinates,
            as well as the intended galactic longitude and current time.
        '''

        # It may be more dangerous to check the connection only once, but
            # heiles is slow and we want to reduce the number of unnecessary calculations
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

    def sweep(self, list_targets):
        start = end = None
        for i in range(len(list_targets)):
            el = list_targets[i][0]
            be = list_targets[i][1]
            alt, az, see = self.find_point_safe(el, be, True)
            if see and start is None:
                start = i
            if see:
                end = i
        if start is None and end is None:
            print('No part of the galactic plane is currently visible.')
        else:
            print('Start:', list_targets[start][0])
            print('End:', list_targets[end][0])

# handy splice:
# list_ell = np.linspace(-10, 250, 261); list_be = np.zeros(261); list_coords = list(zip(list_ell, list_be))
