import ugradio

def g(d):
    """
    Abbreviation function: get
    Acquire data through pico sampler. 1V range,
    divisor=@d.
    """
    return ugradio.pico.capture_data('1V', divisor=10)

def d():
    """
    Abbreviation function: dual get hard-coded for week 2 data collection.
    Acquire two streams of data through pico sampler channels A and B.
    1V range, divisor=2
    """
    return ugradio.pico.capture_data('1V', divisor=2, dual_mode=True)
