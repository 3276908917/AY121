import ugradio
from ugradio import interf_delay
import time
import numpy as np

# always return this to false before committing
dev_mode = False

dc = ugradio.interf_delay.DelayClient()
if dev_mode == False:
    dc.delay_ns(0.)

def jd():
    ''' abbreviation function for current julian date. '''
    return ugradio.timing.julian_date()

def constrain_flip(alt_plan, az_plan):
    '''
    Return adjusted altitude and azimuth so that we take advantage of the
    dish-flip technique for expanding the angular range of the dishes.
    Credit for concept: Rebecca Gore
    '''
    if az_plan <= ugradio.interf.AZ_MIN or \
       az_plan >= ugradio.interf.AZ_MAX:
        # We need to subtract because we are mirroring the
        # altitude across the vertical line 90 degrees
        alt_flip = (180 - alt_plan) % 360
        # We need to add because we are rotating by a
        # a half-circle and not mirroring
        az_flip = (180 + az_plan) % 360
        return alt_flip, az_flip
    return alt_plan, az_plan
    

class Irf:
    def __init__(self, equinox='J2000',
        latitude=ugradio.coord.nch.lat,
        longitude=ugradio.coord.nch.lon,
        altitude=ugradio.coord.nch.alt,
    ):
        '''
        Initialize an interferometery observation object.
        The default parameters are correct for
        New Campbell Hall after the year 2000.
        '''
        self.lat = latitude
        self.lon = longitude
        self.alt = altitude
        self.eq = equinox

        self.ctrl = ugradio.interf.Interferometer()
        self.multi = ugradio.hp_multi.HP_Multimeter()

    def test_system(self):
        '''
        One datum collected by putting the dishes in maintenance mode.
        Credit to Mehdi Ahror.
        '''
        self.ctrl.maintenance()

        data = self.multi.read_voltage()
        print(data)

        self.ctrl.stow()
        return np.savez('test_data', data)

    ### Section: positioning functions. Get altitude and azimuth for...

    def sun(self):
        ''' Return the current position of the Sun. '''
        ra_sun, dec_sun = ugradio.coord.sunpos()
        return ugradio.coord.get_altaz(ra_sun, dec_sun,
            lat=self.lat, lon=self.lon, alt=self.alt, equinox=self.eq)

    # This was an ad-hoc function and should probably be trashed.
    def sun_at(self, j_date):
        ''' Return the position of the Sun at a given julian date j_date. '''
        ra_sun, dec_sun = ugradio.coord.sunpos(j_date)
        return ugradio.coord.get_altaz(ra_sun, dec_sun, j_date,
            lat=self.lat, lon=self.lon, alt=self.alt, equinox=self.eq)
    
    def moon(self):
        ''' Return the current position of the Moon. '''
        ra_moon, dec_moon = ugradio.coord.moonpos(
            lat=self.lat, lon=self.lon, alt=self.alt)
        return ugradio.coord.get_altaz(ra_moon, dec_moon,
            lat=self.lat, lon=self.lon, alt=self.alt, equinox=self.eq)

    # I am not sure about my syntax here, it looks a little fishy.
    # Specifically, I am not sure if the child function will keep track of the variables.
    def star(self, old_ra, old_dec):
        '''
        Return a function which calculates the precessed azimuth and altitude
        for a point source based on its declination and the time of request.
        @old_ra : original right ascension for the epoch
            Convert hour angle directly to an angle via ha(hour, minute, second)
                from rotations.py
            !DO NOT use the intuitive conversion
                right-ascension = LST - hour-angle
        @old_da : original declanation for the epoch
            For convenience of precision, we have
                dec(degree, minute, second) from rotations.py
        '''
        def stargazer():
            prec_ra, prec_dec = ugradio.coord.precess(old_ra, old_dec, equinox=self.eq)
            return ugradio.coord.get_altaz(prec_ra, prec_dec,
                lat=self.lat, lon=self.lon, alt=self.alt, equinox=self.eq)
        return stargazer

    ### end section

    def verify_repos(self, alt_target, az_target, actual):
        '''
        If there is a discrepancy of more than a fifth of a degree,
        between the angles that we desire and the angles that the
        dish reports, we raise an exception.
        '''
        if abs(alt_target - actual['ant_w'][0]) > .2 \
           or abs(az_target - actual['ant_w'][1]) > .2 \
           or abs(alt_target - actual['ant_e'][0]) > .2 \
           or abs(az_target - actual['ant_e'][1]) > .2:
            raise AssertionError('Target is out of range!')

    def reposition(self):
        '''
        The object uses its current coordinate function to re-calculate the
        altitude and azimuth of the target of interest.
        '''
        alt_target, az_target = constrain_flip(self.coord()[0], self.coord()[1])
        self.ctrl.point(alt_target, az_target, wait = True)
        actual = self.ctrl.get_pointing()

        # I still have not confirmed this one
        # verify_repos(alt_target, az_target, actual)
        
        return alt_target, az_target

    # catch keyboard interrupt? low priority issue
    def capture(self, label, total_capture_time = 3960,
                reposition_interval = 60, backup_interval = 600,
                capture_interval = 1, snooze_time=0):
        '''
        @label : string used to prefix each .npz file name
        Credit for original write: Mehdi Arhror
        '''
        # Sleeping risks pipe breaks. The alternative is to start the
        # script regardless of the current time, then to pare down in post.
        time.sleep(snooze_time)
        
        recording_start = last_backup = time.time()
        self.multi.start_recording(capture_interval)
        meta_record = []
        
        while total_capture_time >= time.time() - recording_start :
            # Did we succed in repositioning the dish?
            repos_success = False
            try:
                # Attempt to reposition
                az, alt = self.reposition()
                repos_success = True
            except:
                # If we fail, record the angle we intended
                az, alt = constrain_flip(self.coord()[0], self.coord()[1])

            # meta_record is a parallel array which records
            # the angle (or intended angle) of the dish at a given time
            meta_record.append(
                np.array([az, alt, time.time(), repos_success])
            )
                
            if time.time() - last_backup >= backup_interval:
                data = self.multi.get_recording_data()
                # create a time-stamp for the data
                minutes = str(np.around((time.time() - recording_start) / 60, 2))
                data_name = 'data/' + label + '_' + minutes + '_minutes'
                np.savez(data_name, data=data, stamp=meta_record)
                last_backup = time.time()
    
            time.sleep(reposition_interval)

        self.ctrl.stow()
        np.savez('data/' + label + '_final', data=data, stamp=meta_record)
        self.multi.stop_recording()
        print('No pipe breaks encountered.')
