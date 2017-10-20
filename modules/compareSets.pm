package compareSets;

# sub check_range_overlap($on_end_overlap_y_n, $startA, $stopA, $startB, $stopB)
# sub merge_ovelapping_ranges($array_reference);
# sub check_range_lists_overlap($array_reference_1, $array_reference_2)
# end sub list

########################################################## Universal perl module header ##########################################################

# Load libraries that this module depends on
use warnings;
use strict;
use Sys::Hostname;
use File::Copy;
use Bio::SeqIO;
use Cwd;
use threads;
use diagnostics;

# Set paths to scripts and modules. Setting explicit paths to the scripts and modules in this specific repository (rather than adding paths to @INC, PERLLIB and PATH on
# your local system) avoids the risk of scripts calling the wrong scripts/modules if you have other repositories on your system that happen to have some script- and module names
# in common with this repository.
my $thisfile = (__FILE__);
my $modules = "";
if($thisfile =~ m/^(.+)\//)	{	$modules = $1;	}
my $scripts = $modules;
$scripts =~ s/modules/scripts/;

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

# Get environment information
my ($timestamp, $r, $rscript, $compress, $uncompress) = envir::senseEnv();

# end header

########################################################## Functions ##########################################################

# Checks if two ranges, A and B, overlap each other. Example: A=19-35, B=30-40. A overlaps B and extends it to the left. You can choose whether to define adjacence as overlap (A and B are 'end on end')
# or if you require strict overlap (parameter 'strict_overlap_y_n' is set to 'n' and 'y' respectively)
sub check_range_overlap
	{
	# Set error messages and accept input parameters
        my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
        my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$strict_overlap_y_n, \$startA, \$stopA, \$startB, \$stopB)'\n\nwhere".
	"\t\$strict_overlap_y_n. If this is set to 'n', feature_A is considered to overlap feature_B if it stops directly adjacent to the start of feature_B, and similarly at the end of feature_B.\n".
	"\t\tIf set to 'y' they have to actually overlap.\n".
        "\t\$startA is the start position of feature A\n".
	"\t\$stopA is the stop position of feature A\n".
        "\t\$startB is the start position of feature B\n".
        "\t\$stopB is the stop position of feature B\n\n";

        my @pars = @_ or die $usage;
        foreach my $el (@pars)  {       $el = text::trim($el);  }
	my $strict_overlap = shift @pars or die $usage;
        my $startA = shift @pars or die $usage;	# NB! This will not work if the argument is the number 0 (because it will then be interpreted as false). In that case you need to use 'shift(@pars) or 0', but it may lead to problems....
	my $stopA = shift @pars or die $usage;
	my $startB = shift @pars or die $usage;
	my $stopB = shift @pars or die $usage;

	if($startA > $stopA)	{	die "The stop position of fetaure A is before its start position. Something must be wrong!\n";	}
	if($startB > $stopB)    {       die "The stop position of fetaure B is before its start position. Something must be wrong!\n";  }

	my $startB_adj = $startB;
	my $stopB_adj = $stopB;
	if($strict_overlap eq "n")
		{
		$startB_adj = $startB-1;
		$stopB_adj = $stopB+1;
		}

	my $result="To be determined";

	# If A stops before B starts (A is before B)
	if($stopA < $startB_adj)	{	$result = "A_before_B";	}

	# If A stops no later than B stops	AND	If A starts before B starts (A extends B to the left)
	elsif(($stopA <= $stopB) and ($startA < $startB))	{	$result = "A_extend_B_left";	}

	# If A starts at or after the start of B		AND	If A stops at or before the stop of B (A is contained entirely in B) 
	elsif(($startA >= $startB) and ($stopA <= $stopB))	{	$result = "A_in_B";	}

	# If A starts before B		AND		A stops after B (A extends B in both directions)	
	elsif(($startA < $startB) and ($stopA > $stopB))	{	$result = "A_extend_B_both";	}

	# If A starts before or at the stop of B	AND	If A stops after B (A extends B to the right)
	elsif(($startA <= $stopB_adj) and ($stopA > $stopB))	{	$result = "A_extend_B_right";	}

	# If A starts after B stops (A is after B)
	elsif($startA > $stopB_adj)	{	$result = "A_after_B";	}

	return($result);
	}

# Takes lists
sub merge_ovelapping_ranges
        {
        # Set error messages and accept input parameters
        my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
        my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$array_reference, \$strict_overlap_y_n)'\n\nwhere".
        "\t\$array_reference is a reference to an array holding individual ranges in the format 'start1_stop1, start2_stop2, start3_stop3 etc.\n".
	"\t\$strict_overlap_y_n. If this is set to 'n', feature_A is considered to overlap feature_B if it stops directly adjacent to the start of feature_B, and similarly at the end of feature_B.\n\n";

        my @pars = @_ or die $usage;
        foreach my $el (@pars)  {       $el = text::trim($el);  }
        my $arref = shift @pars or die $usage;      # NB! This will not work if the argument is the number 0 (because it will then be interpreted as false). In that case you need to use 'shift(@pars) or 0', but it may lead to proble$
	my $strict_overlap = shift @pars or die $usage;
	my @arr = @{$arref};	

	# Sort the range array ascending on start and stop
	my @matrix=();
	
	for(my $cc=0; $cc<=$#arr; $cc++)
		{
		my @start_stop = split("_", $arr[$cc]);
		push(@matrix, \@start_stop); 
		}

	@matrix = sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } @matrix;

	my @new_arr=();
	my $comp_start="Aunt_Elfrieda";		# By initially setting this to a string, we will get an error message later if the variable fails to be updated where it should (otherwise we may never know)
	my $comp_stop="Edward_the_Confessor";

	for(my $dd=0; $dd<=$#matrix; $dd++)
		{
		my $start="William";
		my $stop="Harold";
		my $raw_start = $matrix[$dd][0];
		my $raw_stop = $matrix[$dd][1];

		# Just to be sure, check that the start value is indeed smaller than the stop value (and switch them if they are not)
		if($raw_start>$raw_stop)
			{
			$start = $raw_stop;
			$stop = $raw_start;
			}
		else
			{
			$start = $raw_start;
			$stop = $raw_stop;
			}

		# If this is the first line in the matrix, set initial values of $satrt and $stop
		if($dd==0)
			{
			$comp_start = $start;
			$comp_stop = $stop;
			}

		# If it is any other line
		else
			{
			# Check if the current range overlaps the comparison range
			my $overlap=compareSets::check_range_overlap("y", $startA, $stopA, $startB, $stopB);

			}

		}

        return(@new_arr);
	}


# Takes two lists (arrays) each conatining elements that each have a specific location along a discrete scale (the steps ate whole numbers, 1,2,3,4...).
# Determines whiich elements in list A that overlap elements in list B. Creates two outfiles: one containing the elements from list A that overlap
# elements in list B and one with the ones that don't.
sub check_feature_lists_overlap
	{
        # Set error messages and accept input parameters
        my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
        my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$strict_overlap_y_n, \$listA_array_ref, \$listB_array_ref)'\n\nwhere".
        "\t\$strict_overlap_y_n. If this is set to 'n', feature_A is considered to overlap feature_B if it stops directly adjacent to the start of feature_B, and similarly at the end of feature_B.\n".
        "\t\tIf set to 'y' they have to actually overlap.\n".
        "\t\$startA is the start position of feature A\n".
        "\t\$stopA is the stop position of feature A\n".
        "\t\$startB is the start position of feature B\n".
        "\t\$stopB is the stop position of feature B\n\n";

        my @pars = @_ or die $usage;
        foreach my $el (@pars)  {       $el = text::trim($el);  }
        my $strict_overlap = shift @pars or die $usage;
        my $startA = shift @pars or die $usage; # NB! This will not work if the argument is the number 0 (because it will then be interpreted as false). In that case you need to use 'shift(@pars) or 0', but it may lead to problems....
        my $stopA = shift @pars or die $usage;
        my $startB = shift @pars or die $usage;
        my $stopB = shift @pars or die $usage;
	




	}

return(1);

# end functions
