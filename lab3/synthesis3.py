import os

# Grab Lukas' code from labs 1 and 2
# because 2 automatically also grabs 1
os.chdir(r"../lab2")
exec(open('synthesis2.py').read())

# Grab Lukas' code from lab 3

os.chdir(r"../lab3/lukas/src")

exec(open('./plotter.py').read())
exec(open('./start_analysis.py').read())
exec(open('./utils.py').read())

os.chdir(r"../..")
