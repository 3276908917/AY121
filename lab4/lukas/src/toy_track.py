import numpy as np
import ugradio
import ugradio.leusch as leusch

import time
# make sure that you copy rotations.py into the same directory
import rotations

class Plane():
    def __init__(self, latitude=ugradio.leo.lat, longitude=ugradio.leo.lon):
        self.lat = latitude
        self.lon = longitude

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
        import random
        alt_target, az_target, valid = self.find_point_safe(el, be)
        now = ugradio.timing.local_time()
        
        if valid:
            # self.telescope.point
            alt_true = random.getrandbits(8)
            az_true = random.getrandbits(8)
            #self.telescope.get_pointing()

            #self.noise.on()

            #self.lo.set_frequency(634, 'MHz')
            #self.spec.read_spec
            
            #self.lo.set_frequency(635, 'MHz')
            #self.spec.read_spec

            #self.noise.off()
            
            #self.lo.set_frequency(634, 'MHz')
            #self.spec.read_spec
            
            #self.lo.set_frequency(635, 'MHz')
            #self.spec.read_spec
            
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
        # self.spec.check_connected()
        meta_record = []
        
        for coordinate_pair in list_targets:
            el = coordinate_pair[0]
            be = coordinate_pair[1]
            meta_record.append(self.single_measurement(
                el, be, label + '_' + str(el) + '_degrees', N)
            )
            
        np.savez(label + '_stamp', stamp=meta_record)

        print('Ready. If you are done, remember to stow.')

    # I expanded out the default argument
    # to enhance readability.
    def auto_capture(
        self, list_targets, label, N=10,
        sleep_interval = 3600 / 2,
        time_limit = 3600 * 12.5
    ):
        '''
        I want to go to sleep.
        '''
        # self.spec.check_connected()
        meta_record = []
        already = []
        
        start_time = ugradio.timing.unix_time()
        while (ugradio.timing.unix_time() - start_time < time_limit):
            # Do we still have targets to acquire?
            if len(already) < len(list_targets):
                for coordinate_pair in list_targets:
                    if coordinate_pair not in already:
                        el = coordinate_pair[0]
                        be = coordinate_pair[1]
                        attempt = self.single_measurement(
                            el, be, label + '_' + str(el) + '_degrees', N)
                        # Did we actually hit anything?
                        if attempt[4]:
                            meta_record.append(attempt)
                            # This keeps the meta_record from cluttering-up,
                            # because the capture does not save .fits files
                            # anyway if the pointing failed.
                            already.append(coordinate_pair)
                            # Thus, it is important that
                            # list_targets is not a numpy array.
            else:
                break
            # Try again in a little bit, when the firmament has rotated.
            print('\nCycle complete. Sleeping...')
            time.sleep(sleep_interval)
            print('Begin next cycle.\n')
        
        np.savez(label + '_stamp', stamp=meta_record)
        # self.telescope.stow()        

    def sweep(self, list_targets):
        ''' What can I see at the moment? '''
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

# handy splices:
# list_ell = np.linspace(-10, 250, 261); list_be = np.zeros(261); list_coords = list(zip(list_ell, list_be))

# old
# remains = list_coords[:150] + list_coords[190:]

# latest:
# remains = list_coords[:57] + list_coords[129:151] + list_coords[189:]

