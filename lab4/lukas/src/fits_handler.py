from astropy.io import fits

# sandbox / note-taking function
def fits_do(path):
    f = fits.open(path)
    print(f[0].header)
    print(f[1].header)
    data = f[1].data['auto0_real']
    # we sample at 768 MHz but discard 31/32 samples
    x = freq_range(24e6, len(data))
    # 12 MHz Nyquist bandwidth...
        # Do I need to account for 12th window in the
            # x-axis (frequency axis) somehow?

def auto_print(path, index, collector):
    assert collector == 1 or collector == 0, 'Invalid collector index'
    f = fits.open(path)
    auto = f[index].data['auto' + str(collector) + '_real']
    # needs proper frequencies
    plt.plot(auto)

def on_off_print(path_start, path_end, index, collector):
    assert collector == 1 or collector == 0, 'Invalid collector index'

    column = 'auto' + str(collector) + '_real'

    f_on = fits.open(path_start + '_635MHz_' + path_end)
    auto_on = f_on[index].data[column]
    plt.plot(auto_on, label='LO: 1270 MHz')

    f_off = fits.open(path_start + '_634MHz_' + path_end)
    auto_off = f_off[index].data[column]
    plt.plot(auto_off, label='LO: 1268 MHz')

    plt.xlabel('Spectrum Index', fontsize=12)
    plt.ylabel('Power [arbitrary units]', fontsize=12)

    plt.tick_params(axis="x", labelsize=12)
    plt.tick_params(axis="y", labelsize=12)
    plt.legend(bbox_to_anchor=(1, 1))
    
    # needs proper frequencies

def noise_comparison(path_start, index, collector):
    assert collector == 1 or collector == 0, 'Invalid collector index'

    column = 'auto' + str(collector) + '_real'

    f_noise = fits.open(path_start + '_noisy.fits')
    auto_noise = f_noise[index].data[column]
    plt.plot(auto_noise, label='With noise')

    f_q = fits.open(path_start + '_quiet.fits')
    auto_q = f_q[index].data[column]
    plt.plot(auto_q, label='Without Noise')

    plt.xlabel('Spectrum Index', fontsize=12)
    plt.ylabel('Power [arbitrary units]', fontsize=12)

    plt.tick_params(axis="x", labelsize=12)
    plt.tick_params(axis="y", labelsize=12)
    plt.legend(bbox_to_anchor=(1, 1))
    
    # needs proper frequencies

# integration interval comparison
def II_comparison(path, index1, index2, collector):
    assert collector == 1 or collector == 0, 'Invalid collector index'
    f = fits.open(path)

    column = 'auto' + str(collector) + '_real'

    auto_i1 = f[index1].data['auto' + str(collector) + '_real']
    auto_i2 = f[index2].data['auto' + str(collector) + '_real']

    plt.plot(auto_i1, label='Integration interval ' + str(index1 + 1))
    plt.plot(auto_i2, label='Integration interval ' + str(index2 + 1))

    plt.xlabel('Spectrum Index', fontsize=12)
    plt.ylabel('Power [arbitrary units]', fontsize=12)

    plt.tick_params(axis="x", labelsize=12)
    plt.tick_params(axis="y", labelsize=12)
    plt.legend(bbox_to_anchor=(1, 1))
    
    # needs proper frequencies

# absolute data path
adp = 'data/test/'
ending = '_err.npz'

def plane_error_plots(label):
    raw_error = load_saves(adp + label + ending)
    
    alt_on_error = raw_error['on_alt_e']
    az_on_error = raw_error['on_az_e']

    alt_off_error = raw_error['off_alt_e']
    az_off_error = raw_error['off_az_e']
    
    fig, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=True)

    plt.xlabel('Index of Data', fontsize=12)
    fig.text(0, 0.5, r'Angle Discrepancy [degrees]',
             va='center', rotation='vertical', fontsize=12)

    ax1.tick_params(axis="x", labelsize=12)
    ax2.tick_params(axis="x", labelsize=12)
    ax1.tick_params(axis="y", labelsize=12)
    ax2.tick_params(axis="y", labelsize=12)

    ax1.plot(alt_on_error, label='altitude, on')
    ax1.plot(az_on_error, label='azimuth, on')
    ax1.legend(loc='upper right')
    
    ax1.plot(alt_off_error, label='altitude, off')
    ax1.plot(az_off_error, label='azimuth, off')
    ax2.legend(loc='upper right')

# the format is essentially hard-coded
def topo_error_plot(path):
    stamp = load_saves(path)['stamp']

    list_true_alt = stamp[:, 4]
    list_true_az = stamp[:, 5]
    list_target_alt = stamp[:, 2]
    list_target_az = stamp[:, 3]
    list_gal_lat = stamp[:, 0]
    
    fig, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=True)

    plt.xlabel('Galactic Latitude [degrees]', fontsize=12)
    fig.text(0, 0.5, r'Topocentric Angle [degrees]',
             va='center', rotation='vertical', fontsize=12)

    ax1.tick_params(axis="x", labelsize=12)
    ax2.tick_params(axis="x", labelsize=12)
    ax1.tick_params(axis="y", labelsize=12)
    ax2.tick_params(axis="y", labelsize=12)

    ax1.plot(list_gal_lat, list_true_alt, label='True altitude')
    ax1.plot(list_gal_lat, list_target_alt, label='Target altitude')
    ax1.legend(loc='upper right')
    
    ax2.plot(list_gal_lat, list_true_az, label='True azimuth')  
    ax2.plot(list_gal_lat, list_target_az, label='Target azimuth')
    ax2.legend(loc='upper right')
    
