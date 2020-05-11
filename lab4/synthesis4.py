import os

# First, grab Lukas' code from labs 1, 2, and 3 (3 recursively grabs the rest).
os.chdir(r"../lab3")

try:
    exec(open('synthesis3.py').read())
except:
    print('Failed to import lab 1-3 materials.')

# Now grab Lukas' code from lab 4.

os.chdir(r"../lab4/lukas/src")

try:
    exec(open('transform_maps.py').read())
    exec(open('fits_handler.py').read())
    exec(open('velocity.py').read())
    exec(open('mass.py').read())
except:
    print('Failed to import lab 4 materials.')

os.chdir(r"../..")
