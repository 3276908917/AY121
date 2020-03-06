import numpy as np

def fit(L):
    return np.polyfit(null_range(len(L)), L, 1)

def length(L, GHz):
    return fit(L)[0] * GHz * 1e9

# free-space wavelength
def fs(GHz):
    return 3e10 / (GHz * 1e9)

def null_range(L):
    indices = np.zeros(L)
    for i in range(L):
        indices[i] = float(i + 1)
    return indices

# c = \lambda \nu
def guide_velocity(wl, frq):
    return wl * frq

def red_chi(N, M, err, y, yh):
    return 1 / (N - M) / err ** 2 * sum(abs(y - yh) ** 2)     
