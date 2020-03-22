import ugradio

def bounds_pst(stamp_array):
    '''
    Displays (does not return) the start time and stop time for one
    session of data collection based on its associated stamp array.
    '''
    start_ux = stamp_array[:, 2][0]
    start_pst = ugradio.timing.local_time(start_ux)
    print('First datum collected at approximately:\t' + start_pst)
    
    end_ux = stamp_array[:, 2][len(stamp_array) - 1]
    end_pst = ugradio.timing.local_time(end_ux)
    print('Last datum collected at approximately:\t' + end_pst)
