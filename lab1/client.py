# I think I still wanted a routine somewhere that picked out
    # a number of maxima offered by the user

exec(open('synthesis.py').read())                                            

def ACF():
    redos = load_saves('redos.npz')
    redos = load_saves('data/redos.npz')
    r1 = redos['re1']                                                            
    r1c = (r1[100:]).astype(float)
    p1 = power_plot(r1c, 658) 
    i1 = invf(p1)
    plot(i1[0][600:625], np.real(i1[1][600:625]), r'Time ($\mu$s)', 'Voltage (mV)')   
    n_corr = np.correlate(r1c, r1c, mode='full')  
    p_ACF = mini_volt_spec(n_corr, 658) 

# It seems to me like the FT of the ACF looks like the FT of the input signal.
# I am struggling to demonstrate either of the following:
    # Power spectrum looks like FT of ACF
    # iFT of power spectrum looks like ACF

def recreation():
    """
    Keep in mind that this routine was written
    as an emulation of steps in the shell; 
    for example, if you try to refer to 'mix'
    after running this, you will get a not-
    defined error.
    """
    mix = load_saves('data/7_1.npz')
    ld = mix['l_d1'][100:] # axe the bad start part
    PL = power_plot(ld, 435, srate=62.5e6)
    plot(PL[0], PL[1], 'Frequency (MHz)', r'Magnitude-squared Voltage (V$^2$)')

# L Tracker: what can I take hits on?
    # ACF may not be a biggie
    # Most of 7.3

def zeros(P):
    for i in range(0, 32):
        P[i][7990:8008] = 0

def filtering():
    data71 = load_saves('data/7_1.npz')
    ld = data71['l_d1']
    f = ugradio.dft.dft(ld, vsamp=62.5e6)
    # acquire target
    plt.plot(np.real(f[1]))
    # target aquired
    f[1][2400:2600] = f[1][13400:13600] = 0
    # target eliminated
    
def latest_save():
    np.savez('power_pack_7-3.npz', phn = pp_hn, pln = pp_ln, plp = pp_lp, php = pp_hp)

# Low power: sum at 21.4505 MHz
    # difference at .550314 MHz

# High power: sum at 
    # difference at 22.5511 MHz
    # difference at .550314 MHz -> hey, the error is consistent!
