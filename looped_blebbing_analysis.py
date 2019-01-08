# Code to perform looped analysis of membrane blebbing, changing a single parameter to predefined values to examine effect
# Motivation:	it may not always be obvious what the best parameter values are - in this case, 
#				looping over values and looking at the output might be useful
# D.J. Kelly, 2019-01-08, douglas.kelly@riken.jp

# intended operation: have user prepare txt/csv file in format "parameter name (copied from parameters used.json): 
#					[val1, val2, val3]"; when looped analysis is run, prompt for the location of this file and save
#					to output folders named for parameter name and value
