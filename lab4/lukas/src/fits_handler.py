# absolute data path
adp = '../../data/'

def error(label):
    alt_error
    
    fixed_ell = np.array(bound_altaz(ell))
    alts = fixed_ell[:, 0]
    azs = fixed_ell[:, 1]
    
    fig, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=True)

    plt.xlabel('LST [degrees]', fontsize=12)
    fig.text(0, 0.5, r'Angle [degrees]',
             va='center', rotation='vertical', fontsize=12)

    ax1.tick_params(axis="x", labelsize=12)
    ax2.tick_params(axis="x", labelsize=12)
    ax1.tick_params(axis="y", labelsize=12)
    ax2.tick_params(axis="y", labelsize=12)

    ax1.plot(LST, alts, label='altitude')
    ax1.plot(LST, [leusch.ALT_MIN] * len(LST), label='minimum allowed')
    ax1.plot(LST, [leusch.ALT_MAX] * len(LST), label='maximum allowed')
    ax1.legend(loc='upper right')
    
    ax2.plot(LST, azs, label='azimuth')  
    ax2.plot(LST, [leusch.AZ_MIN] * len(LST), label='minimum allowed')
    ax2.plot(LST, [leusch.AZ_MAX] * len(LST), label='maximum allowed')
    ax2.legend(loc='upper right')
