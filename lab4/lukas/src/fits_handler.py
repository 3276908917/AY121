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
