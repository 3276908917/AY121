import ugradio
import time

## aliasing

def jd():
    return ugradio.timing.julian_date()

class Irf:
    def __init__(self, equinox='J2000',
        latitude=ugradio.coord.nch.lat, longitude=ugradio.coord.nch.lon, altitude=ugradio.coord.nch.alt,
    ):
        self.lat = latitude
        self.lon = longitude
        self.alt = altitude
        self.eq = equinox

        self.ctrl = ugradio.interf.Interferometer()
        self.multi = ugradio.hp_multi.HP_Multimeter()

    def test_system(self):
        self.ctrl.maintenance()

        data = self.multi.read_voltage()
        print(data)

        self.ctrl.stow()
        return np.savez('test_data', data)

    # Idea: query user for inputs, then run the data collection routine

    ### Section: positioning functions. Get altitude and azimuth for...

    def sun(self):
        ra_sun, dec_sun = ugradio.coord.sunpos(jd())
        return ugradio.coord.get_altaz(ra_sun, dec_sun, jd(), self.lat, self.lon, self.alt, self.eq)
    
    def moon(self):
        ra_moon, dec_moon = ugradio.coord.moonpos(jd(), self.lat, self.lon, self.alt)
        return ugradio.coord.get_altaz(ra_moon, dec_moon, jd(), self.lat, self.lon, self.alt, self.eq)

    def star(self, old_ra, old_dec):
        '''
        I am not sure about my syntax here, it looks a little fishy.
        Specifically, I am not sure if the child function will keep track of the variables.
        '''
        def stargazer():
            new_ra, new_dec = ugradio.coord.precess(old_ra, old_dec, jd(), self.eq)
            return ugradio.coord.get_altaz(new_ra, new_dec, jd(), self.lat, self.lon, self.alt, self.eq)
        return stargazer

    ### end section
    def point(self):
        alt_target, az_target = self.coord()
        self.ctrl.point(alt_target, az_target, wait = True)
        actual = self.ctrl.get_pointing()
        
        if abs(alt_target - actual['ant_w'][0]) > .1 \
            or abs(az_target - actual['ant_w'][1]) > .1 \
            or abs(alt_target - actual['ant_e'][0]) > .1 \
            or abs(az_target - actual['ant_e'][1]) > .1:
            raise AssertionError('Target is out of range!')
        return True

    def capture(self, label,
                total_capture_time = 3960, reposition_interval = 60,
                backup_interval = 600, capture_interval = 1):
        '''
        ?
        '''
        recording_start = lasst_backup = time.time()

        index = 0
        self.multi.start_recording(capture_interval)

        while total_capture_time >= time.time() - recording_start :
            self.point()
                
            if time.time() - last_backup >= backup_interval:
                data = self.multi.get_recording_data()
                # create a time-stamp for the data
                minutes = str((time.time() - recording_start) / 60)
                data_name = 'data/' + label + '_' + minutes + '_minutes'
                np.savez(data_name, data)
    
            time.sleep(reposition_interval)

        self.ctrl.stow()
        print('No runtime errors encountered.')
