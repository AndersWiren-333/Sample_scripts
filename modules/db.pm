package db;

# select($matrixref, $sortkey, $col, $value);
# union($arref1, $arref2, $num_or_char);
# intersection($arref1, $arref2, $num_or_char);
# end sub list

########################################################## Universal perl module header ##########################################################

# perl_module_update

# Load libraries that this module depends on
use warnings;
use strict;
use Sys::Hostname;
use File::Copy;
use Bio::SeqIO;
use Cwd;
use threads;
#use diagnostics;
use FindBin;
use DBI;
use IO::Compress::Gzip qw(gzip);
use IO::Uncompress::Gunzip qw(gunzip);

# Set paths to scripts and modules. Setting explicit paths to the scripts and modules in this specific repository (rather than adding paths to @INC, PERLLIB and PATH on
# your local system) avoids the risk of scripts calling the wrong scripts/modules if you have other repositories on your system that happen to have some script- and module names
# in common with this repository.
my $thisfile = (__FILE__);
my $modules = "";
if($thisfile =~ m/^(.+)\//)	{	$modules = $1;	}
my $scripts = $modules;
$scripts =~ s/modules/scripts/;
my $maintain = $scripts;
$maintain =~ s/scripts/maintainance/;

# If this script/module is intended to be used outside the folder structure of the parent repository (e.g. a wrapper script to be started from
# another part of your system), set the absolute path to repository scripts and modules (that this cript may depend on) here (and comment out
# the lines for seeting paths above). Otherwise, keep the lines below commented out.
#my $modules = "path/to/modules/folder";
#my $scripts = "path/to/scripts/folder";

# Load home-made modules (aka "libraries", aka "packages") that this module depends on
require "$modules/envir.pm";
require "$modules/dnaTools.pm";
require "$modules/fastaTools.pm";
require "$modules/fastqTools.pm";
require "$modules/matrixTools.pm";
require "$modules/misc.pm";
require "$modules/rtf.pm";
require "$modules/stats.pm";
require "$modules/text.pm";
#require "$modules/compareSets.pm";	# This module is still experimental
require "$modules/fileTools.pm";
require "$modules/combinatorics.pm";
require "$modules/db.pm";
require "$modules/normalise.pm";
require "$modules/listTools.pm";


# Create a timestamp string (can be attached to the name of logfiles, for example
my $timestamp = envir::timestamp();
my $rscript = "Rscript";

# end header

########################################################## Functions ##########################################################


sub select
	{
	# This function takes a (reference to a) matrix ('$matrixref) and selects values from a specified column
	# ($sortkey) for rows where another column ($col) has a specific value ($value). $value can be either
	# numeric or text (set $num_or_char to 'num' or 'char') returns an array with the selected values. If
	# the input matrix has a header, set $header_y_n to "y", otherwise to "n".


	# Set error messages
	my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
	my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$matrixref, \$sortkey, \$col, \$value, \$num_or_char, \$header_y_n)'\n\nwhere".
	"\t\$matrixref is a reference to a matrix holding the values that the selection should be made from\n".
	"\t\$sortkey is the number of the column in the matrix holding the values to be selected\n".
	"\t\$col is the number of the column where values should be \$value for the values in the \$sortkey column to be selected\n".
	"\t\$value is the value that cells in the \$col column should have for the values in the \$sortkey column to be selected\n".
	"\t\$num_or_char is an indicator of whether \$value is numeric (set to 'num') or text (set to 'char')\n".
	"\t\$header_y_n. Set this to 'y' if the input matrix has a header, or to 'n' if it doesn't\n\n";
	
	# Accept input parameters
	my @pars = @_ or die $usage;
	foreach my $el (@pars)  {       $el = text::trim($el);  }
	my $matref = shift @pars or die $usage;	# NB! This will not work if the argument is the number 0 (because it will then be interpreted as false). In that case you need to use 'shift(@pars) or 0', but it may lead to problems....
	my $sortkey = shift @pars or die $usage;
	my $col = shift @pars or die $usage;
	my $value = shift @pars or die $usage;
	my $num_or_char = shift @pars or die $usage;
	my $header_y_n = shift @pars or die $usage;

	if($sortkey eq "n")	{	$sortkey = 0;	}
	if($col eq "n")	{	$col = 0;	}
	if($value eq "n")	{	$value = 0;	}
	
	# Read matrix
	my @matrix=@{$matref};
	if($header_y_n eq "y")	{	my $orig_header_ref = shift(@matrix);	}
	my @outarr=();
	
	# Loop over matrix rows
	for(my $cc=0; $cc<=$#matrix; $cc++)
		{
		if($num_or_char eq "num")	{	if($matrix[$cc][$col] == $value)	{	push(@outarr, $matrix[$cc][$sortkey]);	}	}
		elsif($num_or_char eq "char")	{	if($matrix[$cc][$col] eq $value)	{	push(@outarr, $matrix[$cc][$sortkey]);	}	}
		}
	
	return(@outarr);
	} # end select


sub union
	{
	# Takes (references to) two arrays as input returns the union of their elements, i.e. it combines
	# the two arrays into one and then removes duplicate elements from it.


	# Set error messages
	my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
	my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$arref1, \$arref2, \$num_or_char)'\n\nwhere".
	"\t\$arref1 is a reference to the first array to be used\n".
	"\t\$arref2 is a reference to the second array to be used\n".
	"\t\$num_or_char is an indicator of whether the elements in the two arrays are numeric (set to 'num') or text (set to 'char')\n\n";
	
	# Accept input parameters
	my @pars = @_ or die $usage;
	foreach my $el (@pars)  {       $el = text::trim($el);  }
	my $arref1 = shift @pars or die $usage;	# NB! This will not work if the argument is the number 0 (because it will then be interpreted as false). In that case you need to use 'shift(@pars) or 0', but it may lead to problems....
	my $arref2 = shift @pars or die $usage;
	my $num_or_char = shift @pars or die $usage;
	
	if($num_or_char eq "char")	{	$num_or_char = "alph";	}
	
	my @arr1 = @{$arref1};
	my @arr2 = @{$arref2};
	my @new_arr = (@arr1, @arr2);
	my @unique_array=misc::unique_list(\@new_arr, $num_or_char);

	return(@unique_array);
	} # end union
	

sub intersection
	{
	# Takes (references to) two arrays as input returns the intersection of their elements, i.e. it finds
	# which elements thw two arrays have in common and returns them as an array.


	# Set error messages
	my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
	my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$arref1, \$arref2, \$num_or_char)'\n\nwhere".
	"\t\$arref1 is a reference to the first array to be used\n".
	"\t\$arref2 is a reference to the second array to be used\n".
	"\t\$num_or_char is an indicator of whether the elements in the two arrays are numeric (set to 'num') or text (set to 'char')\n\n";
	
	# Accept input parameters
	my @pars = @_ or die $usage;
	foreach my $el (@pars)  {       $el = text::trim($el);  }
	my $arref1 = shift @pars or die $usage;	# NB! This will not work if the argument is the number 0 (because it will then be interpreted as false). In that case you need to use 'shift(@pars) or 0', but it may lead to problems....
	my $arref2 = shift @pars or die $usage;
	my $num_or_char = shift @pars or die $usage;
	
	if($num_or_char eq "char")	{	$num_or_char = "alph";	}
	my @unique1=misc::unique_list($arref1, $num_or_char);
	my @unique2=misc::unique_list($arref2, $num_or_char);
	
	my @intersection=();
	
	# Loop over rows in @unique1
	ONE: for(my $cc=0; $cc<=$#unique1; $cc++)
		{
		TWO: for(my $dd=0; $dd<=$#unique2; $dd++)
			{
			if($num_or_char eq "num")	{	if($unique1[$cc] == $unique2[$dd])	{	push(@intersection, $unique1[$cc]);	next ONE;	}	}
			elsif($num_or_char eq "alph")	{	if($unique1[$cc] eq $unique2[$dd])	{	push(@intersection, $unique1[$cc]);	next ONE;	}	}
			}
		}

	return(@intersection);
	} # end intersection
	

return(1);

# end functions
