# Description: This is a template for writing new R scripts
# It takes no arguments and does nothing with them
# Usage: Rscript template_script.R
# end description

########################################## Universal R script header ##########################################

# Load libraries (= packages = modules) that this script depends on
library(ggplot2)	# For making graphs. Not part of base distribustion
library(devtools)	# For making packages. Not part of base distribustion
library(knitr)		# For makring R markdown reports. Not part of base distribustion

# Require (load) functions you have written in a separate file
source("my_functions.R")

# Read command line arguments
arguments=commandArgs(trailingOnly=T)
parameter_1=arguments[1]
rest_of_parameters=argu[2:length(arguments)]

# Get date
date=Sys.Date()

# end header

########################################## Define local functions ##########################################

# Declare and define a function (forward declaration - declaring first and defining at the end - is not available in R. I think)
myFunction <- function(argument1=3, argument2=4)
	{
    # This is a template function.
    # It takes two numbers as arguments (default valuea are 3 and 4) and returns the sum of them
    # Usage: myFunction(argument1, argument2)
    result <- argument1+argument2
    return(result)
	}

# end functions

########################################## Processing ##########################################


    

    
    
# end processing