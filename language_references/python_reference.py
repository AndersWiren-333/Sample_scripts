#####     PYTHON QUICK REFERENCE     #####

# Help
# Packages and requiring them
# Scalars
# Arrays
# Arithmetic
# Matrices
# Plots
# Loops
# Conditional statements
# Reading and writing to files
# Functions
# Miscellaneous functions


### HELP ###
help(function)


### PACKAGES AND REQUIRING THEM ###

# Installing packages from command line
python -m pip install --upgrade pip setuptools wheel	# First time, use pip to update/install pip itself, setuptools and wheel
pip install "packagename"		# Istall a new package
pip install --upgrade "packagename"		# Udate an existing package
pip install --user "packagename"		# Istall a new package to the current user's directory/space
pip install -r requirements.txt			# Install all packages listed in "requirements.txt"
# See https://packaging.python.org/tutorials/installing-packages/ for more options

# Load required library/package
import packagename



### SCALARS ###

# Assign a value to a scalar
string_var = "Ale stenar"
num_var = 56

# Splitting a string into an array

### ARRAYS ###

# In python, an array is called a list, unless it has been created with the numpy package, in which case
# it is called an array (or numpy array).

# Assign values to an array (a list)
array = [1,2,3]
array = ["Anders", "Barbro", "Cissi"]
range(1,10,2)		# Makes the array (1,3,5,7,9)

# Get length of an array
len(arrayname)

# Indexing an array
arr = [1,2,3,4,5,6,7,8,9,10]
arr[0]					# 1
arr[-1]					# 10
arr[3]					# 4
arr[2:5]				# 3,4,5
arr[0:10:2]				# Every other element (1,3,5,7,9)
arr[0:len(arr):2]		# Same as above
arr[::2]				# Same as above

# Concatenating arrays
arr3 = arr1+arr2			# [1,2,3,4,5,6,7,8,9,10,1,2,3,4,5,6,7,8,9,10]
arr4 = arr1*3				# [1,2,3,4,5,6,7,8,9,10,1,2,3,4,5,6,7,8,9,10,1,2,3,4,5,6,7,8,9,10]

# Copy an array
ett = [1,2,3,4,5]
tva = ett				# NB! NB! NB! This declaration makes tva a reference to ett! If ett is modified, tva will also be modified! NB! NB! NB!
tva = list(ett)			# This makes an independent copy of the array

# Other array functions
arrayname.append(element)			# Appends one (not more) elements to the ends of an array
del arrayname[2]					# Removes an element from an array
arrayname.reverse()				# Reverses the order of an array
if "Elsa" in namelist				# If the element "Elsa" is present in the array "namelist", do something
if "A" in "Alfabetet"				# If "A" is present in the string "Alfabetet", do something

# Sorting an array
sorted(arrayname)		# Works on both numeric and character arrays


### ARITHMETIC ###
2+1			# 3
3-2			# 1
2*4			# 8
8/2			# 4
2**4 		# 16
3//2		# 1 (get only the integer part of the quotient)

# Get length of a string
len(stringname)

# Concatenating strings
newstring = string1+string2			"Ale stenarAle stenar"
string3 = string1 * 3				"Ale stenarAle stenarAle stenar"

# Indexing a string
string_var[4:5]		# "s"
string_var[:5]		# "Ale s"
string_var[4:]		# "stenar"
string_var[-1]		# "r"
string_var[2:-1]	# "e stena"
string_var[0:10:2]		# Get every other character in the string, "Aesea"
string_var[0:len(string_var):2]		# Same as above
string_var[::2]						# Same as above	


-
### MATRICES ###

# In python, matrices are called 2-dimensional arrays

# Assign values to a matrix
matrix = numpy.array([[1,2,3],[4,5,6],[7,8,9]])

# Read a csv file
import numpy
matrix = numpy.loadtxt(fname="filename.csv", delimiter=",")

# Indexing a matrix
value = matrix[2,4]		# Get the value in the third row and the fifth column (indexing starts at 0)
section = matrix[0:4,0:10]		# Get the first 4(!) rows of the 10(!) first columns ("4" here means up to but excluding 4)
section = matrix[:3,10:]		# Get rows 0, 1 and 2 for all columns from 10 to the last column
array = matrix[2,:]			# Get all columns of the third row and put them in an array
column_as_array = matrix[:,3]		# Gets column 4 and turns it into a 1d array
column_as_column = matrix[:, 3:4]	# Gets column four and preserves it as a column
new_matrix_with_2_cols = numpy.hstack(matrix[:, 3:4], matrix[:, 4:5])		# Stacks two columns horizontally
new_matrix_without_col4 = numpy.delete(matrix, 3, 1)		# Delete one column starting at column 4

# Array/matrix arithmetic
doublematrix = matrix * 2		# Every number in matrix gets multiplied by 2
doublematrix + matrix			# Every number in matrix is added to the corresponding number in doublematrix

# Summary statistics of matrices
import numpy
numpy.mean(matrix)			# Computes the mean of all values in the matrix
numpy.max(matrix)
numpy.min(matrix)
numpy.stdev(matrix)

numpy.mean(matrix, axis=0)		# Computes the average of each column, i.e. for each column, computes the average of all rows for that column
numpy.mean(matrix, axis=1)		# Computes the average of each row, i.e. for each row, computes the average of all columns in that row.

# Get the dimensions of a matrix (or a 1d-array) (rows, cols)
matrix.shape

# Sorting a matrix



### PLOTS ###

# Make plots
import matplotlib.pyplot

heatmap = matplotlib.pyplot.imshow(matrix)		# Make a heat map
lineplot = matplotlib.pyplot.plot(matrix)		# Make a line plot

matplotlib.pyplot.show()						# Show the latest plot
matplotlib.pyplot.show(heatmap)						# Show a specific plot


# Make a figure with subplots
import matplotlib.pyplot

figur = matplotlib.pyplot.figure(figsize=(10.0,3.0))

plot_area1 = figur.add_subplot(1,3,1)			# 1 row, 3 columns, 1st column
plot_area2 = figur.add_subplot(1,3,2)			# 1 row, 3 columns, 2nd column
plot_area3 = figur.add_subplot(1,3,3)			# 1 row, 3 columns, 3rd column

medel = numpy.mean(matrix, axis=0)
maxi = numpy.max(matrix, axis=0)
mini = numpy.min(matrix, axis=0)

plot_area1.set_ylabel('average')
plot_area1.set_ylim(0,10)
plot_area1.plot(medel)

plot_area2.set_ylabel('max')
plot_area2.plot(maxi)

plot_area3.set_ylabel('min')
plot_area3.plot(mini)

figur.tight_layout()

matplotlib.pyplot.show()



### LOOPS ###

# For loop

# While loop

# Foreach loop
for i in list:
	print(i)



### REGULAR EXPRESSIONS ###	

import glob				# This package handles regular expressions in Python

*		# 0 or more of any character
?		# 1 of any character

array = glob.glob("Samples_*_csv")			# Makes an array with filenames (from the current directory) matching the pattern


### CONDITIONAL STATEMENTS ###
number = -5

# if, elif, else
if number > 0:
	print("Greater than 0")
elif number == 0:
	pass
else:
	print("Less than 0")
	
# and, or, not
if (number > 0) and (1 == 1):
	print("-5 is greater than 0 and 1 equals 1")
else:
	print("At least one test is false")
	

if (number > 0) or (1 == 1):
	print("At least one test is true")
else:
	print("Both tests are false")

	
if not 3 == 5:
		print("3 is not 5")
	

# Comparison operators (these work on both strings and numbers)
==	# Equal to
>=	# Greater than or equal to
<=	# Less than or equal to
!=	# Not equal to
	
	
### READING AND WRITING TO FILES ###

# Reading a textfile
# Reading a matrix (table) file

# Writing to a textfile
# Writing to a matrix (table) file
# Writing to a pdf



###  FUNCTIONS  ###

# Definition of user defined function
def function_name(argument1, argument2=0):			# NB! All arguments without a default value must come first in the definition!
	'''Description of what the function does
	Usage: "function_name(argument1, argument2)" where the default for argument2 is 0'''		# This is read by pythons built-in help system. Invoke with "help(function_name)"
	assert precondition, "Precondition not met"		# E.g. check number and type of input arguments
	result = argument1 + argument2
	assert postcondition, "Postcondition not met"		# E.g.check that the return value has the correct type and is within expected range
	return result

# Function call
function_name(1,argument2=2)

# Get help
help(function_name)


# Forward declaration

	
###  MISCELLANEOUS FUNCTIONS ###

# System call

# Time etc.
import time
time.ctime()		# Get the current time and date

# List the contents of a directory

# Assertions
assert 0 < x <=10, "x needs to be larger than 0 and smaller or equal to 10"			# If the condition is not met, the program exits with the message