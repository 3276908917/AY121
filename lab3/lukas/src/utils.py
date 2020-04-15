import matplotlib.pyplot as plt
import ugradio

# Cygnus A "Right Ascension" and Declanation
CAR = ha(19, 59, 28.3566)
CAD = dec(40, 44, 2.096)

def lockout_plotter(stamp):
    fig, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=True)

    times = stamp[:, 6]
    x = times - min(times)

    plt.xlabel('Elapsed Time [s]', fontsize=12)
    fig.text(0, 0.5, r'Angle [degrees]',
         va='center', rotation='vertical', fontsize=12)

    ax1.tick_params(axis="x", labelsize=12)
    ax2.tick_params(axis="x", labelsize=12)
    ax1.tick_params(axis="y", labelsize=12)
    ax2.tick_params(axis="y", labelsize=12)

    alt_want = stamp[:, 0]
    alt_w = stamp[:, 1]
    alt_e = stamp[:, 2]

    ax1.plot(x, alt_want, label='Desired altitude')
    ax1.plot(x, alt_w, label='Actual west dish alt')
    ax1.plot(x, alt_e, label='Actual east dish alt')

    ax1.legend(bbox_to_anchor=(1, 1))
    
    az_want = stamp[:, 3]
    az_w = stamp[:, 4]
    az_e = stamp[:, 5]

    ax2.plot(x, az_want, label='Desired azimuth')
    ax2.plot(x, az_w, label='Actual west dish az')
    ax2.plot(x, az_e, label='Actual east dish az')

    ax2.legend(loc='upper right')

    plt.show()

def bounds_pst(stamp_array, time_index=2):
    '''
    Displays (does not return) the start time and stop time for one
    session of data collection based on its associated stamp array.

    @time_index: which column of the stamp array holds the time?
        This parameter is unfortunately necessary, because the data-
        taking script changed significantly over the course of the lab.
    '''
    start_ux = stamp_array[:, time_index][0]
    start_pst = ugradio.timing.local_time(start_ux)
    print('First datum collected at approximately:\t' + start_pst)
    
    end_ux = stamp_array[:, time_index][len(stamp_array) - 1]
    end_pst = ugradio.timing.local_time(end_ux)
    print('Last datum collected at approximately:\t' + end_pst)

def next_time_boundary(pos_func, start=ugradio.timing.unix_time()):
    alt, az = pos_func(ugradio.timing.julian_date(start))
    add_sec = 60
    # since we want to find the NEXT boundary,
    # we first have to get out of bounds
    while alt > 5 and alt < 175:
        alt, az = pos_func(ugradio.timing.julian_date(start + add_sec))
        add_sec += 60
    if add_sec > 60:
        print('Target sets at unix time', start+add_sec)
        return ugradio.timing.local_time(start + add_sec)
    # now we add time until we are in bounds again
    while alt < 5 or alt > 175:
        alt, az = pos_func(ugradio.timing.julian_date(start + add_sec))
        add_sec += 60
    print('Target rises at unix time', start+add_sec)
    return ugradio.timing.local_time(start + add_sec)

def moon_altaz(j_date):
    ra_moon, dec_moon = ugradio.coord.moonpos(jd=j_date)
    return ugradio.coord.get_altaz(ra_moon, dec_moon, jd=j_date)

def sun_altaz(j_date):
    ra_sun, dec_sun = ugradio.coord.sunpos(jd=j_date)
    return ugradio.coord.get_altaz(ra_sun, dec_sun, jd=j_date)

def define_star(old_ra, old_dec):
    def star_altaz(j_date):
        prec_ra, prec_dec = ugradio.coord.precess(old_ra, old_dec, jd=j_date)
        return ugradio.coord.get_altaz(prec_ra, prec_dec, jd=j_date)
    return stargazer
