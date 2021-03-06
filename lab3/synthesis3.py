import os

# Grab Lukas' code from labs 1 and 2 (2 automatically also grabs 1)
os.chdir(r"../lab2")

try:
    exec(open('synthesis2.py').read())
except:
    print('Failed to import lab 2 materials.')

# Grab Lukas' code from lab 3

os.chdir(r"../lab3/lukas/src")

try:
    exec(open('./analysis3.py').read())
    exec(open('./utils.py').read())
except:
    print('Failed to import lab 3 materials.')

os.chdir(r"../..")
