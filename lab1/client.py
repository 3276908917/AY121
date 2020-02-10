# I think I still wanted a routine somewhere that picked out
    # a number of maxima offered by the user

exec(open('synthesis.py').read())                                            
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


# L Tracker: what can I take the hit on?
    # ACF may not be a biggie
    # Most of 7.3



# Low power: sum at 21.4505 MHz
    # difference at .550314 MHz

# High power: sum at 
    # difference at 22.5511 MHz
    # difference at .550314 MHz -> hey, the error is consistent!
