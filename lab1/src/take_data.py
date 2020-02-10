import ugradio

def g(d):
    """
    Abbreviation function: get
    Acquire data through pico sampler. 1V range,
    divisor=@d.
    """
    return ugradio.pico.capture_data('1V', divisor=10)

def collect_noise():
    """
    Hard coded for final data collection activity of week 1
    """
    noises = []
    for c in range (0, 32):
        noises.append(g(2))
    return noises

def d():
    """
    Abbreviation function: dual get hard-coded for week 2 data collection.
    Acquire two streams of data through pico sampler channels A and B.
    1V range, divisor=2
    """
    return ugradio.pico.capture_data('1V', divisor=2, dual_mode=True)

def power_array(n, norm):
    P = []
    for c in n:
        f = ugradio.dft.dft(c, vsamp=62.5e6)    
        nextP = np.abs(f[1]) ** 2
        x = f[0] / 10 ** 6
        y = normalize(nextP, norm ** 2) / 1000 ** 2
        P.append([x, y])
    return P

# Here is an example of a labeled save: 
# np.savez('bundle1', trial1=t1, trial2=t2, trial3=t3, trial4=t4, trial5=t5, trial6_1=t6_1, trial6_2=t6_2, trial6_3=t6_3, trial7=t7, noises=n) 
