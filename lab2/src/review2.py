def complex_combine(real, imag):
    return real + 1j*imag

# combines power level arrays, chiefly for use with the subsequent
# average and median methods.
def power_lists(P, N):
    stack = [P[i][1] for i in range(1, N)]

def samples_median(power_list):
    return np.median(power_list, axis=0)

def samples_mean(power_list):
    return np.mean(power_list, axis=0)

# "It's okay to degrade the frequency resolution to, say, 1 or 2 kHz'
	# I should probably write a function to handle this.

# Big: calculate reference frame adjustments

def gain(scal, scold):
    return sum(scold) * 300/sum(scal - scold)

# I guess the line frequency is the accepted value for the HI line?
def doppler(nu_0, nu):
    return 3e10 * (nu_0 - nu) / nu_0
