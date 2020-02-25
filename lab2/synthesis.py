import os

os.chdir(r"../lab1")
exec(open('synthesis.py').read())

os.chdir(r"../lab2/src")
exec(open('take_data2.py').read())
exec(open('review2.py').read())
exec(open('analysis2.py').read())

os.chdir(r"../")
