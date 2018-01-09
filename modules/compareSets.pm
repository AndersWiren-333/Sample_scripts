package compareSets;

# check_range_overlap($on_end_overlap_y_n, $startA, $stopA, $startB, $stopB)
# merge_overlap_range($array_reference, $strict_overlap_y_n, $log_y_n);
# check_range_lists_overlap($array_reference_1, $array_reference_2)
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


sub check_range_overlap
	{
	# Checks if two ranges, A and B, overlap each other. Example: A=19-35, B=30-40. A overlaps B and extends it to the left.
	# You can choose whether to define adjacence as overlap (A and B are 'end on end') or if you require strict overlap
	# (parameter 'strict_overlap_y_n' is set to 'n' and 'y' respectively)


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
	} # end check_range_overlap


sub merge_overlap_range
	{
	# Takes (a reference to) and array of ranges in the format "123_456" and checks whether any of the ranges overlap each other.
	# If so it merges those ranges. The updated list of ranges (with some merged) is then returned. Optionally, a log file can be
	# written (option $log_y_n = "y"), e.g. for debugging purposes.


	# Set error messages and accept input parameters
	my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
	my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$array_reference, \$strict_overlap_y_n, \$log_y_n)'\n\nwhere".
	"\t\$array_reference is a reference to an array holding individual ranges in the format 'start1_stop1, start2_stop2, start3_stop3 etc.\n".
	"\t\$strict_overlap_y_n. If this is set to 'n', feature_A is considered to overlap feature_B if it stops directly adjacent to the start of feature_B, and similarly at the end of feature_B.\n".
	"\t\$log_y_n. Set this to 'y' if you want a log to be written (e.g. for debugging) or 'n' if you don't\n\n";

	my @pars = @_ or die $usage;
	foreach my $el (@pars)  {       $el = text::trim($el);  }
	my $arref = shift @pars or die $usage;      # NB! This will not work if the argument is the number 0 (because it will then be interpreted as false). In that case you need to use 'shift(@pars) or 0', but it may lead to proble$
	my $strict_overlap = shift @pars or die $usage;
	my $log_y_n = shift @pars or die $usage;
	my @arr = @{$arref};	

	# Open a log file
	open(my $log, ">>", "log_merge_overlapping_ranges_${timestamp}.txt") or die "Subroutine $subname couldn't create logfile\n";
	
	# Convert @arr to a table with two columns, start and stop
	my @matrix=();
	
	for(my $cc=0; $cc<=$#arr; $cc++)
		{
		my @start_stop = split("_", $arr[$cc]);
		push(@matrix, \@start_stop); 
		}

	# Sort the table ascending on start and stop
	@matrix = sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } @matrix;

	my @new_matrix=();
	my $comp_start="Aunt_Elfrieda";		# If the variable fails to be updated to a number where, we will get an error (= good)
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

		# If this is the first line in the matrix, set initial values of $start and $stop
		if($dd==0)
			{
			$comp_start = $start;
			$comp_stop = $stop;
			}

		# If it is any other line
		else
			{
			# Check if the current range overlaps the comparison range
			my $overlap=compareSets::check_range_overlap("y", $comp_start, $comp_stop, $start, $stop);

			# Evaluate result and merge ranges if appropriate
			
			if($overlap eq "A_before_B")
				{
				# Add A to @new_matrix, set B as new reference range (and loop starts over with next range)
				if($log_y_n eq "y")	{	print $log "${comp_start}_${comp_stop} vs ${start}_${stop}:   A_before_B\n";	}
				push(@new_matrix, [($comp_start, $comp_stop)]);
				$comp_start = $start;
				$comp_stop = $stop;
				if($dd==$#matrix)	{	push(@new_matrix, [($comp_start, $comp_stop)]);	}	# If this is the last line, add it's range to @new_matrix
				}
				
			elsif($overlap eq "A_extend_B_left")
				{
				# Extend stop of A to stop of B, loop starts over with next range
				if($log_y_n eq "y")	{	print $log "${comp_start}_${comp_stop} vs ${start}_${stop}:   A_extend_B_left\n";	}
				$comp_stop = $stop;
				if($dd==$#matrix)	{	push(@new_matrix, [($comp_start, $comp_stop)]);	}	# If this is the last line, add it's range to @new_matrix
				}
				
			elsif($overlap eq "A_in_B")
				{
				# Extend start and stop of A to start and stop of B, loop starts over with next range
				if($log_y_n eq "y")	{	print $log "${comp_start}_${comp_stop} vs ${start}_${stop}:   A_in_B. Could possibly happen\n";	}
				$comp_start = $start;
				$comp_stop = $stop;
				if($dd==$#matrix)	{	push(@new_matrix, [($comp_start, $comp_stop)]);	}	# If this is the last line, add it's range to @new_matrix
				}
				
			elsif($overlap eq "A_extend_B_both")
				{
				# Skip over B and go to next range
				if($log_y_n eq "y")	{	print $log "${comp_start}_${comp_stop} vs ${start}_${stop}:   A_extend_B_both\n";	}
				if($dd==$#matrix)	{	push(@new_matrix, [($comp_start, $comp_stop)]);	}	# If this is the last line, add it's range to @new_matrix
				}
				
			elsif($overlap eq "A_extend_B_right")
				{
				# Set start of A to start of B, loop starts over with next range
				if($log_y_n eq "y")	{	print $log "${comp_start}_${comp_stop} vs ${start}_${stop}:   A_extend_B_right. ************This shouldn't happen*****************\n";	}
				$comp_start = $start;
				if($dd==$#matrix)	{	push(@new_matrix, [($comp_start, $comp_stop)]);	}	# If this is the last line, add it's range to @new_matrix
				}
			
			elsif($overlap eq "A_after_B")
				{
				die "Reference range is after current range. Something must be wrong.\n";
				}
			
			} # Loop over all lines except first one ends

		} # Loop over all lines ends
	
	# Convert @new_matrix to an array of the format "111_222", "333_444"
	my @new_arr = ();
	
	for(my $ee=0; $ee<=$#new_matrix; $ee++)
		{
		my $new_val = $new_matrix[$ee][0] . "_" . $new_matrix[$ee][1];
		push(@new_arr, $new_val);
		}
	
	close($log);
	if($log_y_n eq "n")	{	unlink("log_merge_overlapping_ranges_${timestamp}.txt");		}
	return(@new_arr);
	} # end merge_overlap_range


# sub check_feature_lists_overlap
	# {
	# Takes two lists (arrays) each containing elements that each have a specific location along a discrete
	# scale (the steps are whole numbers, 1,2,3,4...). Determines which elements in list A that overlap elements
	# in list B. Creates two outfiles: one containing the elements from list A that overlap elements in list B
	# and one with the ones that don't. This function is useful for determining e.g. which assembled transcripts
	# overlap previously annotated genes, if their locations on a chromosome/contig is known.


	# # Set error messages and accept input parameters
	# my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
	# my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$strict_overlap_y_n, \$listA_array_ref, \$listB_array_ref)'\n\nwhere".
	# "\t\$strict_overlap_y_n. If this is set to 'n', feature_A is considered to overlap feature_B if it stops directly adjacent to the start of feature_B, and similarly at the end of feature_B.\n".
	# "\t\tIf set to 'y' they have to actually overlap.\n".
	# "\t\$startA is the start position of feature A\n".
	# "\t\$stopA is the stop position of feature A\n".
	# "\t\$startB is the start position of feature B\n".
	# "\t\$stopB is the stop position of feature B\n\n";

	# my @pars = @_ or die $usage;
	# foreach my $el (@pars)  {       $el = text::trim($el);  }
	# my $strict_overlap = shift @pars or die $usage;
	# my $startA = shift @pars or die $usage; # NB! This will not work if the argument is the number 0 (because it will then be interpreted as false). In that case you need to use 'shift(@pars) or 0', but it may lead to problems....
	# my $stopA = shift @pars or die $usage;
	# my $startB = shift @pars or die $usage;
	# my $stopB = shift @pars or die $usage;
	




	# } # end check_feature_lists_overlap


return(1);

# end functions
