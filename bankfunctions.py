#function to read through headers and return first non-header line
def read_through_headers(file,prefix) :
	line = file.readline()
	while line[0:len(prefix)] == prefix :
		line = file.readline()
	return line
	