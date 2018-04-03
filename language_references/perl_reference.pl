#####     PERL QUICK REFERENCE     #####

# Help
# Packages and requiring them
# Print
# Operators
# Function calls
# Local and global variables
# Data types and structures
# Scalars
# Arrays
# Matrices
# Hashes
# Loops
# Conditional statements
# Reading and writing to files
# Writing functions
# System calls
# Regular expressions
# Plots
# Command line arguments and input validation
# Miscellaneous functions



### HELP ###
perldoc keyword


### PACKAGES AND REQUIRING THEM ###

# Load libraries that this script depends on
use warnings;
use strict;
use Getopt::Long qw(GetOptions);	# Load module "Getopt::Long" and make sure you only need to say "GetOptions()" to use that function from it.

# Load user defined modules (aka "libraries", aka "packages") that this script depends on
require "$modules/envir.pm";


# Make a Perl script runnable on unix/linux (start with this line)
#!/usr/bin/perl


### PRINT ###
print "God morgon Världen!";
print("God morgon Världen!")
print $variable;
print "God morgon " . $variable . "! Det ar onsdag idag.";
print("God morgon ", "Världen!");
print "\t\n";		# tab and newline


### OPERATORS ###

=		# Assignment
.		# Concatenation
==, !=, <, >, <=, >=		# Number comparison
eq, ne, lt, gt, le, ge		# String comparison
and, or, not				# Logical
&&, ||, !					# Same as above

# Arithmetic
1+2
5-3
3*5
15/3
2**3		# Exponentiation
2.4e6		# 2.4 * 10**6 (2.4 to the power of 6)
3%2			# Modulo (remainder after division). 3%2=1, 6%2=0, 14%3=2

# Precedence
**
*, /, %
+, -

# Autoincrement
$var++;		# $var = $var + 1;
$var--;		# $var = $var - 1;
$var += 2;		# $var = $var + 2;
$var -= 2;		# $var = $var - 2;
$var *= 2;		# $var = $var * 2;
$var /= 2;		# $var = $var / 2;
$var **= 2;		# $var = $var ** 2;
$var %= 2;		# $var = $var % 2;



### FUNCTION CALLS ###
function(arguments);


### LOCAL AND GLOBAL VARIABLES ###
my $variable;	# Local variable
our $variable;	# Global variable


### DATA TYPES AND STRUCTURES ###

# Data types
Numbers
Strings

# Data structures
Scalars
Arrays
Hashes

### SCALARS ###

# Declare and Assign
my $string_or_number;
my $string = "Faster Tinne";
my $number = 3;

# Index (here used with substring)
my $text = "Studentlagenhet";
my $substring = substr($text, 3, 5);		# "dentl"

# Concatenate
my $text = "Alla " . "fåglar " . "kommit " . "ren";
my $text = "Part 1".
"Part2".
"Part3";

# Iterate over

# Print
print "The value of the \$variable is $variable";

# String functions
$string = lc($string);		# convert to lower case
$string = uc($string);		# convert to upper case
my @array = split("", $string);			# Split into array. ABC becomes ("A", "B", "C")

# Reference
my $scalar_reference = \$scalar;
my $scalar = ${$scalar_reference};
my $scalar = $$scalar_reference;



### ARRAYS ###

# Declare and Assign
my @empty_array;
my @empty_array = ();
my @array1 = (1,2,3,4,5);
my @array2 = ("A", "B", "C", "D", "E");

# Index
my $second_element = $array[1];		# Indexing starts at 0
my $last_element = $array[$#array];		# "$#" = last index of array (NB! This is one less than the number of elements in the array)
my @subset = @array1[1..3];		# 2,3,4

# Concate
my @arr3 = (@arr1, @arr2);

# Iterate over
foreach my $element (@array)	{	print $element;	}
for(my $index=0; $index<=$#array; $index++)	{	print "The value of index $index is $array[$index]";	}

# Print
print("@array");	# With spaces
print(@array);		# Without spaces

# Array functions
push(@array, $new_last_element);		# Adds to end
my $last_element = pop(@array);			# Removes from end
unshift(@array, $new_first_element);
my $first_element = shift(@array);
delete($array[$index]);
my $number_of_elements = scalar(@array);
my @reversed_array = reverse(@array);
my $concatenated_elements_string = join("--", @array);		# $string = "A--B--C--D--E";

# Reference
my @array_reference = \@array;
my @array = @{$array_reference};	# Dereference an array reference
my @array = @$array_reference;

# Sort
my @sorted_array = sort(@array);
my @sorted_numeric_array_ascending = sort { $a <=> $b } @array;
my @sorted_string_array_descending = sort { lc($b) cmp lc($a) } @array;		# "lc" = lowercase, or you may get "Anders Bil Cissi anders bil cissi" as a result



### MATRICES ###

# A matrix is a two-dimensional array. It is constructed by putting array references as elements in an array

# Declare and Assign
my @empty_matrix;
my @empty_matrix = ();
my @matrix=([1,2,3],[4,5,6],[7,8,9]);			# [4,5,6] is a reference to the array (4,5,6)

# Index
my $row_reference = $matrix[1];
my $element = $matrix[$row][$col];

# Concatenate
my @matrix3 = (@matrix1, @matrix2);			# The rows that are @matrix2 will follow the rows that are @matrix1

# Iterate over
foreach my $element (@array)	{	print $element;	}
for(my $index=0; $index<=$#array; $index++)	{	print "The value of index $index is $array[$index]";	}

# Print
print(@array);

# Matrix functions
Since matrices are 2-dimensional arrays, array functions apply to them

# Reference
my $matrix_reference = \@matrix;
my $matrix_reference = ["Anders", "Bertil", "Cissi"];	# Anonymous (unnamed) array
my @matrix = @{$matrix_reference};

# Sort
my @sorted_matrix = sort { lc($a->[3]) cmp lc($b->[3]) } @matrix;		# Sort ascending on (the alphabetic) column 4 of the matrix
my @sorted_matrix = sort { lc($a->[3]) cmp lc($b->[3]) || $b->[1] <=> $a->[1] } @matrix;		# Sort ascending on (the alphabetic) column 4 of the matrix, then descending by
																						# (the numeric) column 2.

																					
# Read csv (table file) into matrix
There is no function for this in the base distribution
																						
																						
### HASHES ###

# Declare and Assign
my %empty_hash;			# Declare empty hash
my %hash1 = ("Country" => "Sweden", "Capital" => "Stockholm", "Currency" => "Krona");
my %hash2 = ("Key 1" => 2.14, "Key 2" => 3.5, "Key 3" => 11);
my %hash1 = (									# Also a valid syntax
			"Country" => "Sweden",
			"Capital" => "Stockholm",
			"Currency" => "Krona",
			);
my %hash1 = ("Country", "Sweden", "Capital", "Stockholm", "Currency", "Krona");		# Also a valid syntax
my @array = ("Country", "Sweden", "Capital", "Stockholm", "Currency", "Krona");
my %hash = @array;		# Also a valid syntax

# Index
my $land = $hash1{"Country"};		# "Sweden"
$hash1{"Currency"} = "Euro";		# Change the value of an existing key
$hash1{"Head of state"} = "Carl Gustaf";		# Add a new key-value pair

# Concatenate

# Iterate over
for my $key in keys(%hash)	{	print "The value of the key $key is $hash{$key}";	}

# Print
There is no function for this in the base distribution of perl

# Hash functions
my $key_list = keys(%hash);
my $value_list = values(%hash);
delete($hash{key});

# Reference
my $hash_reference = \%hash;
my $hash_reference = {"Alva" => 23, "Kurt" => 67};		# Anonymous (unnamed) hash
my %hash = %{$hash_reference};
my %hash = %$hash_reference;

# Sort
A hash is an unordered list by design, and I haven't found a way to save a sorted hash,
but a hash can be printed in sort order like this:

foreach my $key (sort { lc($a) cmp lc($b) } keys(%hash))	{	print "$key\t$hash{$key}\n";	}			# Sort alphabetically on key
foreach my $key (sort { $hash{$a} <=> $hash{$b} } keys(%hash))	{	print "$key\t$hash{$key}\n";	}		# Sort numerically on value
foreach my $key (sort { $hash{$a} <=> $hash{$b} || lc($a) cmp lc($b) } keys(%hash))	{	print "$key\t$hash{$key}\n";	}	# Sort numerically on value,
																															# then alphabetically on key



### LOOPS ###

# For loop
for(my $i=0; $i<=5; $i++)	{	print "$i";	}

# Foreach loop
foreach my $el (@array)	{	print "$el";	}

# While loop
my $index=1;
while($index<=10)	{	print "$index"; $index++;	}


###   CONDITIONAL STATEMENTS   ###
if(condition)	{	do this;	}
elseif(condition2)	{	do this;	}
else	{	do default option;	}
unless(condition)	{	do something;	}
print ($verbose ? "$verbose is true" : "$verbose is false");	# The ternary operator. If the variable is true, do the first option, else the other



### READING AND WRITING TO FILES ###

# Open a file handle
open(my $in, "<", "infile.txt") or die "Script $0 couldn't open infile $infile\n";		# Read from file
open(my $out, ">", $outname) or die "Script $0 couldn't create outfile $outname\n";		# For writing
open(my $out, ">>", $outname) or die "Script $0 couldn't create outfile $outname\n";		# For appending to end of open file

# Close a filehandle
close($in);

# Read a directory
my $dir = shift // ".";		# If no directory name given at the command line, use the current directory ("//" means "defined or" - if shift doesn't return a defined value, use ".")
opendir(my $dh, $dir) or die "Script $0 couldn't open directory $dir for reading\n";
my @files = readdir($dh);
# or
my $files = grep { $_ ne "." and $_ ne ".." } readdir($dh);

# or
while(my $file = readdir($dh))	{	print "$file\n";	}
closedir($dh);

# Other file functions
unlink $file		# Remove
use File::Copy; move($old_name, $new_name);			# Rename


### WRITING FUNCTIONS ###
sub my_function;		# Forward declaration
sub my_function			# Definition of function
	{
	my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
	my @pars = @_ or die $usage;	# Rather than @ARGV, which holds command line arguments for a script
	my $result = 1+2;
	return($result);
	}
my_function();		# Function call (= use the function)


### SYSTEM CALLS ###



### REGULAR EXPRESSIONS ###



### PLOTS ###

Use R for this


### COMMAND LINE ARGUMENTS AND INPUT VALIDATION ###

# Simple and without dependencies
my @pars = @ARGV or die $usage;
foreach my $el (@pars)  {       $el = text::trim($el);  }		# Make sure there are no whitespaces at the end or beginning of command line arguments
my $infile = shift(@pars) or die $usage;
my $namefile = shift(@pars) or die $usage;
my @rest_of_pars = @pars or die $usage;

# Using Getopt::Long (base distribution)
# THIS PACKAGE CAN'T HANDLE ARRAY REFERENCES AS ARGUMENTS AND CAN'T PERFORM INPUT VALIDATION - That is why I don't use it.
use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_getopt);		# Enable the use of short names (-v -l etc) and bundling of them (-vl). This can't be used with repeat modifiers (e.g. "{2}")
my $person = "Bishop";			# "Bishop" is the default value
my $place = "Ely";
my $adjective = "a good man";
my $verbose;					# This variable doesn't have a default value
my $repeat=1;
my $charm=1.0;
my @infiles=();
GetOptions(
	"person|p=s" => \$person,			# option "person" (which can be specified as "--person" or "-p") - if it is given - must have a value and the value must be s string
	"place|L:s" => \$place,				# The option - if given - doesn't have to be given a value (if the value is left out, the default will be set to "")
	"adjective|a=s" => \$adjective,
	"verbose|v" => \$verbose,
	"repeat|r=i" => \$repeat,			# The value of the option must be an integer (if ":" instead of "=", default is set to 0)
	"charmfactor|c=f" =>	\$charm,	# The value must be a floating point number (= a decimal number) (if ":" instead of "=", default is set to 0)
	"infiles=s{1,}" => \@infiles,			# Take one or more infiles  ({2,4} = between two and four values, {,} = zero or more values, {2} = exactly two values). Can't be used with bundling.
	) or die "Usage: perl $0 [--person STRING] [--place STRING] [--adjective STRING] [--verbose] [--repeat INT] [--charm FLOAT] [--infiles STRING1 STRING2...]\n";

if($verbose)	{	print "The $person of $place is $adjective.\n" x $repeat;	}
else	{	print "$person $place $adjective" x $repeat;	}

# Data::FormValidator is another package, that I haven't tied yet. Read more about it at
# http://search.cpan.org/~dfarrell/Data-FormValidator-4.88/README.pod



### MISCELLANEOUS FUNCTIONS ###

# Get system information
my $script = $0;		# Name of this (running) script
my $os = $^O;			# Name of operating system
use Sys::Hostname;		# Name of computer that script is running on
my $computer = hostname;
my $dec=`perl -MPOSIX=locale_h -e "print localeconv()->{decimal_point}`;	# Get the decimal separator for the current system

# Compress and uncompress files
use IO::Compress::Gzip qw(gzip);
use IO::Uncompress::Gunzip qw(gunzip);
gzip "example.txt" => "example.txt.gz";
unlink("example.txt");						# The uncompressed file has to be deleted manually
gunzip "example.txt.gz" => "example.txt";
unlink("example.txt");						# The compressed file has to be deleted manually

# Free up memory when no longer needing a variable/data structure
undef @matrix;
@matrix = ();		# This may not free all the memory, use the syntax above instead

# File tests
if(-e $filename)	{	print "File $filename exists\n";	}
if(-f $filename)	{	print "File $filename is a file\n";	}
if(-d $filename)	{	print "File $filename is a directory\n";	}
if(-l $filename)	{	print "File $filename is a symbolic link\n";	}
if(-z $filename)	{	print "File $filename is empty\n";	}

if(-e $filename && -f _ and -z _)	{	print "File $filename exists and is a file and is empty\n";	}
# or
if(-e -f -z $filename )	{	print "File $filename exists and is a file and is empty\n";	}


# Stop running the script (put this at the end)
exit;

