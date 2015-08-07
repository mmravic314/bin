# Functions for random formatting
## Marco Mravic August 2015 DeGrado Lab UCSF

# Take string of a common delimited list and return a new line delimited list
def comma_2_newLine_list(string):

	new_str = ''
	for i in string.split(','):
		new_str += i.strip() + '\n'

	return new_str

