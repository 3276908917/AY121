import math

G = 6.673e-11 # Gravitational constant in mks
Rb = 4.4e7# estimation of Milky Way supermassive black hole radius, km

def macerate(doppler):
    '''
    Return a radius-dependent mass cumulative distribution function
    based on a single Doppler pair (ell, Doppler_velocity).
    Units are (degrees, km / s).
    '''
    v = velocify(doppler)
    return lambda r: r * v(r) ** 2 / G

def velocify(doppler):
    '''
    Return a radius-dependent velocity function based on
    a single Doppler pair (ell, Doppler_velocity).
    Units are (degrees, km / s).
    '''
    sin_ell = np.sin(np.radians(doppler[0]))
    v_dopp = doppler[1]
    return lambda r: r / R0 * (v_dopp / sin_ell - V0)

# I do not want to learn how to use function objects.
def integrate(func, start, stop, delta):
    '''
    Return the Riemann sum for a function with label @func
    that accepts only one parameter.
    We sum over a linear range of inputs beginning with
    @start and ending at @stop.
    Instead of an infinitesimal, we use @delta to
    provide a finite width to each Riemann rectangle.
    The procedure is specifically a midpoint Riemann sum,
        to avoid tricky things like ill-defined boundary behavior.
    
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
