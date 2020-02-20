import os

os.chdir(r"../lab1")
exec(open('synthesis.py').read())

os.chdir(r"../lab2/src")
exec(open('gestalt.py').read())
exec(open('analysis2.py').read())

os.chdir(r"../")
