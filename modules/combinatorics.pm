package combinatorics;

# sub pairwise_combinations($array_reference, $num_or_char)
# sub remove_self_combinations($matrix_reference, $num_or_char)
# sub remove_redundant_combinations($matrix_reference, $num_or_char)
# sub value_in_array($value, $arrayref)
# sub replace_in_array($old_value, $new_value, $arrayref)
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

# Create a timestamp string (can be attached to the name of logfiles, for example
my $timestamp = envir::timestamp();
my $rscript = "Rscript";

# end header

########################################################## Functions ##########################################################

# Takes a list (an array reference) of elements, numeric or as strings, finds the set of unique values (e.g. 3,5,3,7,4  becomes 3,4,7), and returns
# as a matrix a list of all possible pairwise combinations of these values. Each row in the output matrix represents a combination (e.g. 3,3 or 3,7
# or 7,4). Values can be either numeric or text (set $num_or_char to 'num' or 'char' at command line). Optionally, self-combinations (e.g. 3,3, "A","A")
# can be removed (set option $remove_self_y_n to "y", otherwise "n") from the new matrix before it is returned, as can duplicate combinations (i.e. when 3,4 and 4,3
# is considered the same, one of them is removed - set option $remove_duplicates_y_n to "y", otherwise "n").
sub pairwise_combinations
	{
	# Set error messages and accept input parameters
	my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
	my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$array_ref, \$num_or_char, \$remove_self_y_n, \$remove_duplicates_y_n)'\n\nwhere".
	"\t\$array_ref is a reference to the array holding the values to be combined\n".
	"\t\$num_or_char. Set this to 'num' if the values to be combined are numeric, or to 'char' if they are text\n".
	"\t\$remove_self_y_n. Set this to 'y' to remove self-combinations (e.g. 3,3 or 'A','A') from output matrix, otherwise set it to 'n'\n".
	"\t\$remove_duplicates_y_n. Set this to 'y' to remove duplicates (e.g. 4,3 is a duplicate of 3,4) from output matrix, otherwise set it to 'n'\n\n";

	my @pars = @_ or die $usage;
	foreach my $el (@pars)  {       $el = text::trim($el);  }
	my $arref = shift @pars or die $usage;	# NB! This will not work if the argument is the number 0 (because it will then be interpreted as false). In that case you need to use 'shift(@pars) or 0', but it may lead to problems....
	my $num_or_char = shift @pars or die $usage;
	if(($num_or_char ne "num") and ($num_or_char ne "char"))        {       die "The \$num_or_char argument for sub $subname must be either 'num' or 'char'. Try again!\n";      }
	my $rem_self = shift @pars or die $usage;
	my $rem_dupl = shift @pars or die $usage;
	
	my @pairwise_matrix=();

	# Get list of unique values
	my @raw_unique_array=();
	my @unique_array=();
	
	if($num_or_char eq "num")
		{
		@raw_unique_array=misc::unique_list($arref, "num");
		@unique_array = sort { $a <=> $b } @raw_unique_array;
		}
	elsif($num_or_char eq "char")
		{
		@raw_unique_array=misc::unique_list($arref, "alph");
		@unique_array = sort { $a cmp $b } @raw_unique_array;
		}
	
	undef @raw_unique_array;
	
	# Loop over elements in unique list
	for(my $c=0; $c<=$#unique_array; $c++)
		{
		# Loop over elements in unique list again
		for(my $d=0; $d<=$#unique_array; $d++)
			{
			push(@pairwise_matrix, [($unique_array[$c], $unique_array[$d])]);
			}
		}
	undef @unique_array;
		
	my @selfless_matrix=();
	if($rem_self eq "y")	{	@selfless_matrix=combinatorics::remove_self_combinations(\@pairwise_matrix, $num_or_char);	}
	else	{	@selfless_matrix = @pairwise_matrix;	}
	undef @pairwise_matrix;
	
	my @unique_matrix=();
	if($rem_dupl eq "y")	{	@unique_matrix=combinatorics::remove_redundant_combinations(\@selfless_matrix, $num_or_char);	}
	else	{	@unique_matrix = @selfless_matrix;	}
	undef @selfless_matrix;
	
	return(@unique_matrix);
	}
	
# Takes a table (as a matrix reference) of pairwise combinations of values (two columns) and removes combinations
# of equal values (e.g. 3,3 and "A,A"). Returns a new two column table (as a matrix). Values can be either numeric
# or textual (set $num_or_char to 'num' or 'char' at command line).
sub remove_self_combinations
	{
	# Set error messages
	my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
	my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$matrix_ref, \$num_or_char)'\n\nwhere".
	"\t\$matrix_ref is a reference to the matrix holding the combinations of values to be evalueted\n".
	"\t\$num_or_char. Set this to 'num' if the values to be combined are numeric, or to 'char' if they are text\n\n";

	# Accept input parameters
	my @pars = @_ or die $usage;
	foreach my $el (@pars)  {       $el = text::trim($el);  }
	my $matref = shift @pars or die $usage;	# NB! This will not work if the argument is the number 0 (because it will then be interpreted as false). In that case you need to use 'shift(@pars) or 0', but it may lead to problems....
	my $num_or_char = shift @pars or die $usage;
	if(($num_or_char ne "num") and ($num_or_char ne "char"))        {       die "The \$num_or_char argument for sub $subname must be either 'num' or 'char'. Try again!\n";      }
	
	my @pairwise_matrix = @{$matref};
	
	# Remove self-self combinations (e.g. 3,3 or Alfa,Alfa) from matrix
	my @raw_outmatrix=();

	# Loop over @pairwise_matrix
	for(my $c=0; $c<=$#pairwise_matrix; $c++)
		{
		if($num_or_char eq "num")
			{
			unless($pairwise_matrix[$c][0] == $pairwise_matrix[$c][1])	{	push(@raw_outmatrix, $pairwise_matrix[$c]);	}
			}
		elsif($num_or_char eq "char")
			{
			unless($pairwise_matrix[$c][0] eq $pairwise_matrix[$c][1])	{	push(@raw_outmatrix, $pairwise_matrix[$c]);	}
			}
		}

	my @outmatrix=();
	if($num_or_char eq "num")	{	@outmatrix = sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } @raw_outmatrix;	}
	elsif($num_or_char eq "char")	{	@outmatrix = sort { $a->[0] cmp $b->[0] || $a->[1] cmp $b->[1] } @raw_outmatrix;	}
	
	return(@outmatrix);
	}	

# Takes a table (as a matrix reference) of pairwise combinations of values (two columns) and removes combinations
# that are duplicates of other combinations in the sense that which value comes first in a pair doesn't matter,
# for example, 3,4 and 4,3 are duplicates, as is "A,B" and "B,A". Returns a new two column table (as a matrix).
# Values can be either numeric or textual (set $num_or_char to 'num' or 'char' at command line).
sub remove_redundant_combinations
	{
	# Set error messages
	my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
	my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$matrix_ref, \$num_or_char)'\n\nwhere".
	"\t\$matrix_ref is a reference to the matrix holding the combinations of values to be evalueted\n".
	"\t\$num_or_char. Set this to 'num' if the values to be combined are numeric, or to 'char' if they are text\n\n";

	# Accept input parameters
	my @pars = @_ or die $usage;
	foreach my $el (@pars)  {       $el = text::trim($el);  }
	my $matref = shift @pars or die $usage;	# NB! This will not work if the argument is the number 0 (because it will then be interpreted as false). In that case you need to use 'shift(@pars) or 0', but it may lead to problems....
	my $num_or_char = shift @pars or die $usage;
	if(($num_or_char ne "num") and ($num_or_char ne "char"))        {       die "The \$num_or_char argument for sub $subname must be either 'num' or 'char'. Try again!\n";      }
	
	my @pairwise_matrix = @{$matref};
	
	# Remove duplicate combinations (e.g. 3,4 vs 4.3 or Alfa,Beta vs Beta,Alfa) from matrix
	
	# Loop over @pairwise_matrix and sort elements within rows
	for(my $c=0; $c<=$#pairwise_matrix; $c++)
		{
		my @arr = @{$pairwise_matrix[$c]};
		my @new_arr=();
		if($num_or_char eq "num")	{	@new_arr = sort { $a <=> $b } @arr;	}
		elsif($num_or_char eq "char")	{	@new_arr = sort { $a cmp $b } @arr;	}
		
		$pairwise_matrix[$c] = [@new_arr];
		}

	# Sort sorted matrix on columns
	my @smatrix=();
	if($num_or_char eq "num")	{	@smatrix = sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } @pairwise_matrix;	}
	elsif($num_or_char eq "char")	{	@smatrix = sort { $a->[0] cmp $b->[0] || $a->[1] cmp $b->[1] } @pairwise_matrix;	}

	# Loop over rows of sorted matrix and remove duplicates
	my @outmatrix=();
	my @comparr = @{$smatrix[0]};
	push(@outmatrix, $smatrix[0]);
	
	for(my $d=1; $d<=$#smatrix; $d++)
		{
		my @row = @{$smatrix[$d]};
		unless(@row ~~ @comparr)
			{
			push(@outmatrix, $smatrix[$d]);
			@comparr = @{$smatrix[$d]};
			}
		}
	
	return(@outmatrix);
	}

# Checks whether a value (supplied as a scalar) is present in array (supplied as an array reference) or not. Returns “y” if it is and “n” if it isn’t.
# Set argument $num_or_char to ‘num’ if the value and array are numeric or ‘to char’ if they are text.
sub value_in_array
	{
	# Set error messages
	my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
	my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$value, \$arrayref, \$num_or_char)'\n\nwhere".
	"\t\$value is the value whose presence in the array should be checked\n".
	"\t\$arrayref is a reference to the array in which the presence of \$value should be checked\n".
	"\t\$num_or_char  Set this to 'num' if the \$value and \@array are numeric, or 'char' if they are text\n\n";
	
	# Accept input parameters
	my @pars = @_ or die $usage;
	foreach my $el (@pars)  {       $el = text::trim($el);  }
	my $value = shift @pars or die $usage;	# NB! This will not work if the argument is the number 0 (because it will then be interpreted as false). In that case you need to use 'shift(@pars) or 0', but it may lead to problems....
	my $arref = shift @pars or die $usage;
	my $num_or_char = shift @pars or die $usage;

	my @arr = @{$arref};
	
	# Loop over elements in @arr and check if they corresponds to $value
	my $y_or_n = "n";

	LINES: for(my $c=0; $c<=$#arr; $c++)
		{
		if($num_or_char eq "num")	{	if($arr[$c] == $value)	{	$y_or_n = "y"; last LINES;	}	}
		elsif($num_or_char eq "char")	{	if($arr[$c] eq $value)	{	$y_or_n = "y"; last LINES;	}	}
		}

	return($y_or_n);
	}
	

# Checks whether all the elements in array1 (supplied as an array reference) are present in array2 (also supplied as an array reference) or not.
# Returns “y” if they are and “n” if they aren’t. Set argument $num_or_char to ‘num’ if the arrays are numeric or ‘to char’ if they are text.
sub array_in_array
	{
	# Set error messages
	my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
	my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$arrayref1, \$arrayref2, \$num_or_char)'\n\nwhere".
	"\t\$arrayref1 is a reference to the array the presence of whose values should be checked in the other array\n".
	"\t\$arrayref2 is a reference to the array in which the presence ofthe values of array1 should be checked\n".
	"\t\$num_or_char  Set this to 'num' if the they arrays are numeric, or 'char' if they are text\n\n";
	
	# Accept input parameters
	my @pars = @_ or die $usage;
	foreach my $el (@pars)  {       $el = text::trim($el);  }
	my $arref1 = shift @pars or die $usage;	# NB! This will not work if the argument is the number 0 (because it will then be interpreted as false). In that case you need to use 'shift(@pars) or 0', but it may lead to problems....
	my $arref2 = shift @pars or die $usage;
	my $num_or_char = shift @pars or die $usage;

	my @arr1 = @{$arref1};
	my @arr2 = @{$arref2};
	
	# Loop over elements in @arr1 and check if they are present in @arr2
	my $y_or_n_global = "y";

	ARR1_VAL: for(my $c=0; $c<=$#arr1; $c++)
		{
		my $value = $arr1[$c];
		my $y_or_n=combinatorics::value_in_array($value, $arref2, $num_or_char);
		if($y_or_n eq "n")	{	$y_or_n_global = "n"; last ARR1_VAL;	}
		}

	return($y_or_n_global);
	}

return(1);

# end functions
