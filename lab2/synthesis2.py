import os

os.chdir(r"./../lab1")
exec(open('synthesis1.py').read())

os.chdir(r"./../lab2/src")
exec(open('./data_taking.py').read())
exec(open('./data_appraisal.py').read())
exec(open('./analysis2.py').read())
exec(open('./rotations.py').read())
exec(open('./unpack.py').read())
exec(open('./chi.py').read())

os.chdir(r"./../")
