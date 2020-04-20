    def collect(self, el, be, full_prefix, N):
        
        alt, az = self.find_point(el, be)
        self.telescope.point(alt, az)

        alt_err, az_err = self.pos_error(alt, az)
        alt_real, az_real = self.telescope.get_pointing()
        
        list_alt_err.append(alt_err)
        list_az_err.append(az_err)

        self.spec.check_connected()
        ra, dec = gal_to_eq(el, be)
        self.spec.read_spec(full_prefix + '.fits', N, (ra, dec), 'eq')

        return list_alt_err, list_az_err

    # needs some major revision before we move on to large consecutive pointings
        # ability to point at multiple coordinates
        # taking with noise on and with noise off
    def take_data(self, el, be, label, N=10):
        '''
        Collect @N spectra
        by observing galactic coordinates (@el, @be)
            (currently handles only degrees)
        and save the data in two files, each with @label in the name
        the .fits file stores the spectra
        the .npz file stores the two angle-errors for the pointings.
        '''

        self.noise.off()
        
        self.lo.set_frequency(634, 'MHz')
        on_alt_err, on_az_err = self.collect(el, be, label + '_634MHz_quiet', N)
        self.lo.set_frequency(635, 'MHz')
        off_alt_err, off_az_err = self.collect(el, be, label + '_635MHz_quiet', N)

        self.noise.on()
        
        self.lo.set_frequency(634, 'MHz')
        on_alt_err, on_az_err = self.collect(el, be, label + '_634MHz_noisy', N)
        self.lo.set_frequency(635, 'MHz')
        off_alt_err, off_az_err = self.collect(el, be, label + '_635MHz_noisy', N)

        np.savez(label + '_err',
                 on_alt_e=on_alt_err, on_az_e=on_az_err,
                 off_alt_e=off_alt_err, off_az_e=off_az_err)

        self.noise.off()
        self.telescope.stow()
