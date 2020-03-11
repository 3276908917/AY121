# not tested at all

import os

os.chdir(r"../../lab1")
exec(open('synthesis.py').read())

os.chdir(r"../lab2")
exec(open('synthesis.py').read())

os.chdir(r"src")
# Code files go here

os.chdir(r"./../")
