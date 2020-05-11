import math

Rb = 4.4e7# estimation of Milky Way supermassive black hole radius, km

def velocify(doppler):
    sin_ell = np.sin(np.radians(doppler[0]))
    v_dopp = doppler[1]
    return lambda r: r / R0 * (v_dopp / sin_ell - V0)

# I do not want to learn how to use function objects.
def integrate(func, start, stop, delta):
    '''
    This only works for analytical, continuous functions.
    '''
    num_quanta = math.ceil((stop - start) / delta)
    x_space = np.linspace(start, stop, num_quanta)

    riemann = 0
    for i in range(len(x_space) - 1):
        mid_x = 0.5 * (x_space[i] + x_space[i + 1])
        mid_y = func(mid_x)
        riemann += mid_y * delta

    return riemann
