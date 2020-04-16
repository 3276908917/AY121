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