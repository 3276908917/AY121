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

    def test_system(self, sampling_frequency):
        self.controller.maintenance()

        data = self.multi.read_voltage()
        print (data)

        self.ctrl.stow()
        return np.savez('test_data', data)

# We can worry about user utilities later. Right now we just need some data
   # def user_record(self):
   #     tt = input('For how many hours would you like to record?')
   #     print(tt == '')

    def fast_record(self, coord_func, total_time=3960, reposition_interval=60, backup_interval=600, dt=1):
        '''
        @coords is one parameter, because it is much easier to pass in the tuple that
            the abbreviation functions returns than to calculate them oneself
        '''
        self.user_record()

    ### Section: positioning functions. Get altitude and azimuth for...

    def sun(self):
        ra_sun, dec_sun = ugradio.coord.sunpos(jd())
        return ugradio.coord.get_altaz(ra_sun, dec_sun, jd(), self.lat, self.lon, self.alt. self.eq)

    def moon(self):
        ra_moon, dec_moon = ugradio.coord.moonpos(jd(), self.lat, self.lon, self.alt)
        return ugradio.coord.get_altaz(ra_moon, dec_moon, jd(), self.lat, self.lon, self.alt. self.eq)

    def star(self, old_ra, old_dec):
        '''
        I am not sure about my syntax here, it looks a little fishy.
        Specifically, I am not sure if the child function will keep track of the variables.
        '''
        def stargazer(self):
            new_ra, new_dec = ugradio.coord.precess(old_ra, old_dec, jd(), self.equinox)
            ugradio.coord.get_altaz(new_ra, new_dec, jd(), self.lat, self.lon, self.alt, self.eq)
        return stargazer

    ### end section

    def capture(self, label, coord,
                total_capture_time = 3960, reposition_interval = 60,
                backup_interval = 600, capture_interval = 1):
        '''
        Capture data of the sun every @dt seconds
        '''
        recording_start = lasst_backup = time.time()

        index = 0
        self.multi.start_recording(dt)

        while total_time >= time.time() - recording_start :
            alt_target, az_target = coord()
            self.ctrl.point(alt_target, az_target, wait=True)

            if time.time() - last_backup >= backup_interval:
                data = self.multi.get_recording_data()
                # create a time-stamp for the data
                minutes = str((time.time() - recording_start) / 60)
                data_name = 'data/' + label + '_' + minutes + '_minutes'
                np.savez(data_name, data)
    
            time.sleep(reposition_interval)

        self.ctrl.stow()
        print('No runtime errors encountered.')

    def Record_star(self, old_ra, old_dec, recording_time, total_observation_time):
        count = 0

        ifm = self.initialize_control()
        hpm = self.initialize_voltage()

        ifm.point(Alt, azimuth)

        intitial_time = time.time()

        while time.time() <= total_observation_time + intitial_time:
            print(ifm.get_pointing())

            while time.time() <= initial_time + recording_time:
                time.sleep(time_delay)
                hpm.start_recording(recording_time)
                data = hpm.get_recording_data()
                data_name = 'data_moon_'+ str(count)
                np.save(data_name, data)
                hpm.end_recording()
                count += 1

        ifm.stow()

        return count
