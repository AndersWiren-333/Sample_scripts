package listTools;

# array_equality($arrayref1, $arrayref2)
# get_closest_values($value)
# get_array_value_index($array_reference, $num_or_char, $value)
# get_closest_indices($array_ref, $value)
# value_in_array($value, $array_ref, $num_or_char)
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
# another part of your system), set the absolute path to repository scripts and modules (that this script may depend on) here (and comment out
# the lines for setting paths above). Otherwise, keep the lines below commented out.
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


sub array_equality
	{# Checks whether two arrays (specified as array references) are equal, and returns “y” if they are
	# and “n” if they aren’t.
	
	
	# Set error messages
	my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
	my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$array_reference1, \$array_reference2)'\n\nwhere".
	"\t\$array_reference1 is a reference to the first array to be compared\n".
	"\t\$array_reference2 is a reference to the second array to be compared\n\n";
	
	# Accept input parameters
	my @pars = @_ or die $usage;
	foreach my $el (@pars)  {       $el = text::trim($el);  }
	my $arref1 = shift @pars or die $usage;	
	my $arref2 = shift @pars or die $usage;

	# Processing
	my @arr1 = @{$arref1};
	my @arr2 = @{$arref2};
	my $equal="y";
	
	if(scalar(@arr1) != scalar(@arr2))	{	$equal = "n";	}
	
	else
		{
		# Loop over elements in @arr1 and check if the corresponding element in @arr2 is equal
		for(my $cc=0; $cc<=$#arr1; $cc++)
			{
			# If the two corresponding values are not of the same data type, conclude that the arrays are not equal
			my $type1=misc::type($arr1[$cc]);
			my $type2=misc::type($arr2[$cc]);
			if($type1 ne $type2)	{	$equal = "n";	}
			
			# If they are of the same type, check if the are the same
			else
				{
				if($type1 eq "num")	{	if($arr1[$cc] != $arr2[$cc])	{	$equal = "n";	}	}
				elsif($type1 eq "char")	{	if($arr1[$cc] ne $arr2[$cc])	{	$equal = "n";	}	}
				}
			}
		}

	return($equal);
	} # end array_equality


sub get_closest_values
	{
	# This subroutine checks whether a specific numeric value is present in an array, and if it isn’t it returns the next lowest value
	# and the next highest value in the same array. If the value is present, both the next lowest and highest values are returned as
	# the value itself. If the value is lower than the lowest value in the array, the return values are “” and the smallest value in
	# the array. If the value is larger than the largest value in the array, the return values are the largest value in the array and “”.
	
	
	# Set error messages
	my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
	my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$array_ref, \$value)'\n\nwhere".
	"\t\$array_ref is a reference to the array to be searched in\n".
	"\t\$value is the value to be searched for\n\n";
	
	# Accept input parameters
	my @pars = @_ or die $usage;
	foreach my $el (@pars)  {       $el = text::trim($el);  }	
	my $arref = shift @pars or die $usage;
	my $value = shift @pars or 0;
	my @arr = @{$arref};

	# Processing
	@arr = sort { $a <=> $b } @arr;
	my $low="";
	my $high="";
	
	if($value<$arr[0])
		{
		$low = "";
		$high = $arr[0];
		print "subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}) finds that the".
		" input value is smaller than the smallest value in the array it is being compared to\nThis may cause problems for some applications\n";
		}
	elsif($value>$arr[$#arr])
		{
		$low = $arr[$#arr];
		$high = "";
		print "subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}) finds that the".
		" input value is larger than the largest value in the array it is being compared to\nThis may cause problems for some applications\n";
		}
	
	else
		{
		for(my $cc=0; $cc<=$#arr; $cc++)
			{
			if($arr[$cc]<$value)	{	next;	}
			elsif($arr[$cc]==$value)	{	$low=$value; $high=$value; last;	}
			elsif($arr[$cc]>$value)
				{
				$low = $arr[${cc}-1];
				$high = $arr[$cc];
				last;
				}
			}
		}
	return($low, $high);
	} # end get_closest_values


sub get_array_value_index
	{
	# Checks whether a value is present in an array (specified as an array reference) and if it is,
	# returns the position (index) of that value in the array (0 is the first element in an array).
	# If the value isn’t present, the subrroutine returns “not_available”.

	
	# Set error messages
	my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
	my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$array_reference, \$num_or_char, \$value)'\n\nwhere".
	"\t\$array_reference is a reference to the array to be searched\n".
	"\t\$num_or_char  Set to 'num' if \$value is numeric or to 'char' if it is text\n".
	"\t\$value is the value whose index to search for in the array\n\n";
	
	# Accept input parameters
	my @pars = @_ or die $usage;
	foreach my $el (@pars)  {       $el = text::trim($el);  }
	my $arref = shift @pars or die $usage;
	my $num_or_char = shift @pars or die $usage;
	my $value = shift @pars or 0;
	my @arr = @{$arref};

	# Loop over elements in @arr
	my $index="not_available";
	for(my $cc=0; $cc<=$#arr; $cc++)
		{
		if($num_or_char eq "num")	{	if($arr[$cc] == $value)	{	$index = $cc;	}	}
		elsif($num_or_char eq "char")	{	if($arr[$cc] eq $value)	{	$index = $cc;	}	}
		}
	
	return($index);
	} # end get_array_value_index


sub get_closest_indices
	{
	# This subroutine checks whether a specific numeric value is present in an array, and if it isn’t it returns the
	# index of the next lowest value and the index of the next highest value in the same array. If the value is present,
	# both the next lowest and highest indices are returned as the index of the value itself. If the value is lower than
	# the lowest value in the array, the return values are “” and the index of the smallest value in the array. If the
	# value is larger than the largest value in the array, the return values are the index of the largest value in the
	# array and “”.
	
	
	# Set error messages
	my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
	my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$array_ref, \$value)'\n\nwhere".
	"\t\$array_ref is a reference to the array to be searched in\n".
	"\t\$value is the value whise closest indices are to be searched for\n\n";
	
	# Accept input parameters
	my @pars = @_ or die $usage;
	foreach my $el (@pars)  {       $el = text::trim($el);  }	
	my $arref = shift @pars or die $usage;
	my $value = shift @pars or 0;
	my @arr = @{$arref};

	# Processing
	my @rank_matrix=();
	my @orig_order=(1..scalar($#arr));
	@rank_matrix=matrixTools::add_cols_to_matrix(\@rank_matrix, \@orig_order, \@arr);

	@rank_matrix = sort { $a->[1] <=> $b->[1] } @rank_matrix;
	my $low_rank="";
	my $high_rank="";
	
	# If $value is smaller than the smallest value in the array
	if($value<$rank_matrix[0][1])
		{
		$low_rank = "";
		$high_rank = ($rank_matrix[0][0])-1;
		print "subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}) finds that the".
		" input value is smaller than the smallest value in the array it is being compared to\nThis may cause problems for some applications\n";
		}
		
	# If $value is larger than the largest value in the array
	elsif($value>$rank_matrix[$#rank_matrix][1])
		{
		$low_rank = ($rank_matrix[$#rank_matrix][0]);
		$high_rank = "";
		print "subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}) finds that the".
		" input value is larger than the largest value in the array it is being compared to\nThis may cause problems for some applications\n";
		}
	
	# If value is within the range of the array
	else
		{
		for(my $cc=0; $cc<=$#rank_matrix; $cc++)
			{
			if($rank_matrix[$cc][1]<$value)	{	next;	}
			elsif($rank_matrix[$cc][1]==$value)
				{
				$low_rank = ($rank_matrix[$cc][0])-1;
				$high_rank = ($rank_matrix[$cc][0])-1;
				last;
				}
			elsif($rank_matrix[$cc][1]>$value)
				{
				$low_rank = ($rank_matrix[$cc][0])-2;
				$high_rank = ($rank_matrix[$cc][0])-1;
				last;
				}
			}
		}
	return($low_rank, $high_rank);
	} # end get_closest_indices


sub value_in_array
	{
	# Checks whether a value (supplied as a scalar) is present in array (supplied as an array reference)
	# or not. Returns 'TRUE' if it is and 'FALSE' if it isn’t. Set argument $num_or_char to ‘num’ if the value
	# and array are numeric or to 'char' if they are text.


	# Set error messages
	my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
	my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$value, \$arrayref, \$num_or_char)'\n\nwhere".
	"\t\$value is the value whose presence in the array should be checked\n".
	"\t\$arrayref is a reference to the array in which the presence of \$value should be checked\n".
	"\t\$num_or_char  Set this to 'num' if the \$value and \@array are numeric, or 'char' if they are text\n\n";
	
	# Accept input parameters
	my @pars = @_ or die $usage;
	foreach my $el (@pars)  {       $el = text::trim($el);  }
	my $value = shift @pars or die $usage;
	my $arref = shift @pars or die $usage;
	my $num_or_char = shift @pars or die $usage;

	my @arr = @{$arref};
	
	# Loop over elements in @arr and check if they corresponds to $value
	if($value eq "zero")	{	$value=0;	}
	my $tf = '';

	LINES: for(my $c=0; $c<=$#arr; $c++)
		{
		if($num_or_char eq "num")	{	if($arr[$c] == $value)	{	$tf = 'TRUE'; last LINES;	}	}
		elsif($num_or_char eq "char")	{	if($arr[$c] eq $value)	{	$tf = 'TRUE'; last LINES;	}	}
		}

	return($tf);
	} # end value_in_array	

	
return(1);

# end functions
