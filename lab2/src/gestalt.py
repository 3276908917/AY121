def d2():
    """
    Abbreviation function: dual get hard-coded for week 1 data collection.
    Acquire two streams of data through pico sampler channels A and B.
    10V range, divisor=2
    """
    return ugradio.pico.capture_data('50mV', divisor=1, dual_mode=True)

def complex_combine(real, imag=None, defaults=True):
	if defaults:
		r = real[100:16000]
		i = real[16100:32000]
		return r + 1j*i
	return real + 1j*imag
