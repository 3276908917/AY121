# TODO:
	# upload pictures from phone so that
	# lab partners can perform visual
	# verifications for themselves


# if you haven't already written these, write them 
	# Fourier plotting method,
	# which labels the real and imaginary
	# and includes a legend on the plot
    # no_title option for plots, so that you can create
    # pictures suitable for the report (where titles
        # will not be tolerated, only captions)

# QOL
	# a global label called 'temp'
	# which all functions edit before
	# returning, so that
	# when I forget to save a result,
	# I can go back and write temp to something else

	# expand the docstrings

# Collaboration
	# put all bundles and 'readme.md' describing each item in the repo
	# now you have just one nexus to put on the Google Doc
	# this will also allow you to maintain a changelog,
	# to consolidate all your memoranda and postscripts

# this file in particular
# to-do: why am I executing these files instead of importing them?
    # because I do not want to have to put the module name in front
    # of the function every time that I call it / refer to it
# to-do, long term: only import what you need for any given session

"""
import take_data
import review
import analysis
"""

exec(open('src/take_data.py').read())
exec(open('src/review.py').read())
exec(open('src/analysis.py').read())

# to-do: implement the temp variable?
