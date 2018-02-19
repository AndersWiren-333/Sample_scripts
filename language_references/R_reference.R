# This is a reference collection of R commands

# Packages
install.packages("analysisTools")	# Install
remove.packages("analysisTools", "path/to/package")		# Uninstall
library("analysisTools")				# Load
detach(package:analysisTools)			# Unload
source("my_functions.R")		# Load functions from an R file (without having to turn them into a package)

# Read command line arguments
arguments=commandArgs(trailingOnly=T)
parameter_1=arguments[1]
rest_of_parameters=argu[2:length(arguments)]

# Get date
date=Sys.Date()

# Declare and define a function (forward declaration - declaring first and defining at the end - is not available in R. I think)
functionName <- function(argument1=3, argument2=4)
	{
    # This is a template function
    # It takes two numbers as arguments (default values are 3 and 4) and returns the sum of them
    # Usage: functionName(argument1, argument2)
    result <- argument1+argument2
    return(result)
	}

# Declare variables and assign values to them
vector <- c(1,2,3,4,5)
sequence <- seq(5)

# For loop
for(i in 1:5)
	{
    print(i)
	}
	
# While loop

# For each loop

# Conditional statement
if((1>0) & (3==3))
  {
  print("Hej")
  } else {          # It is important that both braces and the "else" are on the same line
  print("Då")
  }

# Join strings
new_string <- paste("First", "Second", sep=" ")	

# Set up a pdf to print something to
pdf("something.pdf", height=4, width=10)
dev.off()       # This closes the pdf
    
# Plots
plot(vector)
boxplot(vector)
stripchart(vector)
hist(vector)
   
# Run an R script from the command line
Rscript script.R

# Run an R script from the command line if Rscript doesn't work
R --slave --no-restore --file=script.R --args
    
# Get memory size of object
object.size(object)

# Force R to release memory it is no longer using (R doesn't always do this automatically although you have deleted objects)
gc()

# List all objects in current workspace
ls()

# Delete all objects in current work space
rm(list = ls())

# Read csv file
data <- read.csv(file="infile.csv", header = FALSE)

# Indexing dataframes (matrices)
value <- $data[row, column]
row <- data[rownnumber,]
column <- data[,columnnumber]
column <- data$columnname
    
# end processing