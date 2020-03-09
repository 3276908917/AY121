import ugradio 
import numpy as np
import matplotlib.pyplot as plt

def load_saves(filename):
	a = np.load(filename)
	return dict(zip(("{}".format(k) for k in a), (a[k] for k in a)))
