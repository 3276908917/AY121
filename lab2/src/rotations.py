import numpy as np
import ugradio

eq2g = np.array([[-.054876, -.873437, -.483835], [.494109, -.444830, .746982], [-.867666, -.198076, .455984]]

eq_to_ha(LST):
    return np.array([np.cos(LST), np.sin(LST), 0], [-np.sin(LST), np.cos(LST), 0], [0, 0, 1])

topocentric(az, alt):
    return np.array([np.cos(az) * np.cos(alt), np.sin(az) * np.cos(alt), np.sin(alt)])

# Not yet implemented
gal_to_topo():
