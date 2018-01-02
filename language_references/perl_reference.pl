# Make a oerl script runnable on unix/linux (start with this line)
#!/usr/bin/perl

# Load libraries that this script depends on
use warnings;
use strict;

# Load home-made modules (aka "libraries", aka "packages") that this script depends on
require "$modules/envir.pm";


# Declare a scalar variable and assign a value to it
my $string = "Faster Tinne";
my $number = 3;

# Declare an array amd assign values to it
my @tex_tarray = ("Anders", "Barbro", "Cecilia");
my @numeric_array = (1,2,3);

# Concatenate strings
my $text = "Alla " . "fåglar " . "kommit " . "ren";

# Print a text
print "The value of this $variable will be printed in this message since it is enclosed in double quotes. This is called interpolation";

# Concatenate long texts (that span several lines in the editor)
my $text = "Part 1".
"Part2".
"Part3";

# Accept command line arguments
my @pars = @ARGV or die $usage;
foreach my $el (@pars)  {       $el = text::trim($el);  }		# Make sure there are no whitespaces at the end or beginning of command line arguments
my $infile = shift(@pars) or die $usage;
my $namefile = shift(@pars) or die $usage;
my @rest_of_pars = @pars or die $usage;

# For loop
for(my $i=0; $i<=5; $i++)	{	print "$i";	}

# Foreach loop
foreach my $el (@array)	{	print "$el";	}

# While loop
my $index=1;
while($index<=10)	{	print "$index"; $index++;	}

# Writing functions (= subroutines)
sub my_function;		# Forward declaration
sub my_function			# Definition of function
	{
	my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
	my @pars = @_ or die $usage;	# Rather than @ARGV, which holds command line arguments for a script
	my $result = 1+2;
	return($result);
	}
my_function();		# Function call (= use the function)


# Open a file handle
open(my $in, "<", "infile.txt") or die "Script $0 couldn't open infile $infile\n";		# Read from file
open(my $out, ">", $outname) or die "Script $0 couldn't create outfile $outname\n";		# For writing
open(my $out, ">>", $outname) or die "Script $0 couldn't create outfile $outname\n";		# For appending to end of open file

# Close a filehandle
close($in);

# Get system information
my $script = $0;		# Name of this (running) script
my $os = $^O;			# Name of operating system
use Sys::Hostname;		# Name of computer that script is running on
my $computer = hostname;

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

# Read a directory
my $dir = shift // ".";		# If no directory name given at the command line, use the current directory ("//" means "defined or" - if shift doesn't return a defined value, use ".")
opendir(my $dh, $dir) or die "Script $0 couldn't open directory $dir for reading\n";
my @files = readdir($dh);
# or
my $files = grep { $_ ne "." and $_ ne ".." } readdir($dh);

# or
while(my $file = readdir($dh))	{	print "$file\n";	}
closedir($dh);

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

