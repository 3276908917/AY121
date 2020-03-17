import ugradio
import time
import numpy as np

## aliasing

def jd():
    return ugradio.timing.julian_date()

def constrain_flip(alt_plan, az_plan):
    '''
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
        self.lat = latitude
        self.lon = longitude
        self.alt = altitude
        self.eq = equinox

        self.ctrl = ugradio.interf.Interferometer()
        self.multi = ugradio.hp_multi.HP_Multimeter()

    def test_system(self):
        '''
        Credit to Mehdi Ahror
        '''
        self.ctrl.maintenance()

        data = self.multi.read_voltage()
        print(data)

        self.ctrl.stow()
        return np.savez('test_data', data)

    ### Section: positioning functions. Get altitude and azimuth for...

    def sun(self):
        '''
        Calculate the position of the sun at the time of request.
        '''
        ra_sun, dec_sun = ugradio.coord.sunpos(jd())
        return ugradio.coord.get_altaz(ra_sun, dec_sun, jd(), self.lat, self.lon, self.alt, self.eq)

    # This was an ad-hoc function and should probably be trashed.
    def sun_at(self, jd):
        '''
        Calculate the position of the sun at a given julian date.
        '''
        ra_sun, dec_sun = ugradio.coord.sunpos(jd)
        return ugradio.coord.get_altaz(ra_sun, dec_sun, jd, self.lat, self.lon, self.alt, self.eq)
    
    def moon(self):
        '''
        Calculate the position of the moon at the time of request.
        '''
        ra_moon, dec_moon = ugradio.coord.moonpos(jd(), self.lat, self.lon, self.alt)
        return ugradio.coord.get_altaz(ra_moon, dec_moon, jd(), self.lat, self.lon, self.alt, self.eq)

    # I am not sure about my syntax here, it looks a little fishy.
    # Specifically, I am not sure if the child function will keep track of the variables.
    def star(self, old_ra, old_dec):
        '''
        Return a function which calculates the precessed azimuth and altitude
        for a point source based on its declination and the time of request.
        @old_ra : original right ascension for the epoch
        @old_da : original declanation for the epoch
        '''
        def stargazer():
            new_ra, new_dec = ugradio.coord.precess(old_ra, old_dec, jd(), self.eq)
            return ugradio.coord.get_altaz(new_ra, new_dec, jd(), self.lat, self.lon, self.alt, self.eq)
        return stargazer

    ### end section

    def reposition(self):
        '''
        The object uses its current coordinate function to re-calculate the
        altitude and azimuth of the target of interest.
        '''
        alt_target, az_target = constrain_flip(self.coord()[0], self.coord()[1])
        self.ctrl.point(alt_target, az_target, wait = True)
        actual = self.ctrl.get_pointing()

        # I think this may have been broken because of competition for the dishes
        #if abs(alt_target - actual['ant_w'][0]) > .2 \
        #    or abs(az_target - actual['ant_w'][1]) > .2 \
        #    or abs(alt_target - actual['ant_e'][0]) > .2 \
        #    or abs(az_target - actual['ant_e'][1]) > .2:
        #    raise AssertionError('Target is out of range!')
        return alt_target, az_target

    def capture(self, label, total_capture_time = 3960,
                reposition_interval = 60, backup_interval = 600,
                capture_interval = 1, snooze_time=0):
        '''
        @label : string used to prefix each .npz file name

        Credit for original write: Mehdi Arhror
        '''
        time.sleep(snooze_time)
        
        recording_start = last_backup = time.time()
        self.multi.start_recording(capture_interval)
        coord_record = []
        
        while total_capture_time >= time.time() - recording_start :
            try:
                az, alt = self.reposition()
            except:
                az, alt = constrain_flip(self.coord()[0], self.coord()[1])
            coord_record.append(np.array([az, alt, time.time()]))
                
            if time.time() - last_backup >= backup_interval:
                data = self.multi.get_recording_data()
                # create a time-stamp for the data
                minutes = str(np.around((time.time() - recording_start) / 60, 2))
                data_name = 'data/' + label + '_' + minutes + '_minutes'
                np.savez(data_name, data=data)
                last_backup = time.time()
    
            time.sleep(reposition_interval)

        self.ctrl.stow()
        np.savez('data/' + label + '_final', data=data)
        self.multi.stop_recording()
        print('No runtime errors encountered.')
