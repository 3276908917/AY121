import numpy as np
import scipy as sp
import ugradio
import matplotlib.pyplot as plt
import time

## end section

class Irf:
    def __init__(self,
        julian_day=ugradio.timing.julian_date(), latitude=ugradio.coord.nch.lat,
        longitude=ugradio.coord.nch.lon, altitude=ugradio.coord.nch.alt, equinox='J2000'
    ):
        self.jd = julian_day
        self.lat = latitude
        self.lon = longitude
        self.alt = altitude
        self.eq = equinox

        self.ctrl = ugradio.interf.Interferometer()
        self.multi = ugradio.hp_multi.HP_Multimeter()

    def fourier_frequency(self, signal, sampling_frequency):
        '''
        TODO: this function does NOT plot the real
        and imaginary components of the spectrum,
        so the current output is a cast!
        '''
        signal_fft = np.fft.fft(signal)
        time_step = 1 / sampling_frequency
        freqs = np.fft.fftfreq(signal_fft.size, time_step)
        index = np.argsort(freqs)

        plt.plot(freqs[index] / 1e6, ps[index])
        plt.xlabel('Frequency (MHz)')
        plt.ylabel('Amplitude (V)')
        plt.title('Fourier transform')

        plt.legend()
        plt.show()

        return freqs[index], signal_fft[index]

    def power_spectrum(self, signal):
        fourier_transform = np.fft.fft(signal)
        return np.abs(fourier_transform) ** 2

    def spectrum_frequency(self, signal, sampling_frequency):
        ps = power_spectrum(signal)
        time_step = 1 / sampling_frequency
        freqs = np.fft.fftfreq(signal.size, time_step)
        index = np.argsort(freqs)

        #plot the frequency power spectrum
        plt.plot(freqs[index] / 1e6, ps[index])
        plt.xlabel('Frequency (MHz)')
        plt.ylabel(r'Amplitude ($V^2$)')
        plt.title('Power spectrum')
        # plt.yscale('log')

        plt.legend()
        plt.show()

        return freqs[index], ps[index]

    def test_system(self, sampling_frequency):
        self.controller.maintenance()

        data = self.multimeter.read_voltage()
        print (data)

        self.controller.stow()
        return np.savez('test_data', data)


    def user_record(self):
        tt = input('For how many hours would you like to record?')
        print(tt == '')

    def fast_record(self, total_time=3960, reposition_interval=60, backup_interval=600, dt=1, coords)
        '''
        @coords is one parameter, because it is much easier to pass in the tuple that
            the abbreviation functions returns than to calculate them oneself
        '''

    ## Section: abbreviation functions. Get altitude and azimuth for...

    def sun(self)
        ra_sun, dec_sun = ugradio.coord.sunpos(self.julian_day)
        return ugradio.coord.get_altaz(ra_sun, dec_sun, self.jd, self.lat, self.lon, self.alt. self.eq)

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

    # Ra, Dec = ugradio.coord.moonpos(self.julian_day, self.latitude, self.longitude, self.altitude)

    def Record_star(self, old_ra, old_dec, recording_time, total_observation_time):
        count = 0

        ifm = self.initialize_control()
        hpm = self.initialize_voltage()

        new_ra, new_dec = ugradio.coord.precess(old_ra, old_dec, self.julian_day, self.equinox)

        Alt, azimuth = ugradio.coord.get_altaz(new_ra, new_dec, self.julian_day, self.latitude, self.longitude, self.altitude, self.equinox)

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
