import ugradio

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


    def user_record(self):
        tt = input('For how many hours would you like to record?')
        print(tt == '')

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

    def Record_sun(self, total_time=3960, reposition_interval=60, backup_interval=600, dt=1):
        '''
        Capture data of the sun every @dt seconds
        '''

        ifm = self.initialize_control()
        hpm = self.initialize_voltage()
        initial_time = time.time()

        index = 0
        hpm.start_recording(recording_time)

        while self.Total_recording_time >= time.time() - initial_time :
            ra_target, dec_target = ugradio.coord.sunpos(self.julian_day)

            alt_target, az_target = ugradio.coord.get_altaz(
                ra_target, dec_target, self.julian_day, self.latitude,
                self.longitude, self.altitude, self.equinox
            )

            ifm.point(Alt, azimuth, wait=True)

            if self.data_saved_time * index <= time.time() - initial_time:
                data = hpm.get_recording_data()
                data_name = 'data/data_sun_' + str(index + 1)
                np.savez(data_name, data)
                index += 1

                time_1 = time.time()

            while self.Update_position_time >= time.time() - time_1: # it will read data until it's time to switch position ( every 1 minute is better i guess)
                time.sleep(1)
                             # record 1 data every 'recording_time' seconds


    #			the data will be saved every time the telescope change position (1 minute) not sure if it's convinient, probably 10 mins is better, so i concidered saving the data every 10th time that we change our position.

        ifm.stow()

        return index

    # Ra, Dec = 

    def Record_star(self, old_ra, old_dec, recording_time, total_observation_time):
        count = 0

        ifm = self.initialize_control()
        hpm = self.initialize_voltage()

        

        Alt, azimuth = 

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
