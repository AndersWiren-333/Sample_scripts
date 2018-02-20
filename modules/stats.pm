package stats;

# rand_between($min, $max, $decimals)
# rand_array($n, $min, $max, $decimals)
# mean_above($arrayref, $abund_limit, $samp_limit)
# mean($arrayref)
# num_above_cutoff($arrayref, $cutoff)
# median($arrayref)
# min($arrayref)
# max($arrayref)
# variance($arrayref)
# stdev($arrayref)
# percentile($arrayref, $k)
# ofc($array_ref1, $array_ref2)
# list_overlap($array_ref1, $array_ref2);
# test_diff_expression_range($expression_matrix.csv, groups.txt, $outname)
# geometric_mean($array_ref)
# three_point_expr_pattern()
# assign_ranks($array_ref, $num_or_char)
# covariance($array_ref1, $array_ref2)
# corr_coeff($array_ref1, $array_ref2)
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


sub rand_between
	{
	# Generates a random number between given values and with a given number of decimals
	# Parameters: $min, $max, $decimals
	
	# Set error messages and accept input parameters
	my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
	my $usage = "Usage  error for sub '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$min, \$max, \$decimals);'\n\nwhere".
	"\t\$min is the lower limit of the number to be returned\n".
	"\t\$max is the upper limit of the number to be returned\n".
	"\t\$decimals is the number of desired decimals of the returned number\n\n";
	my $min = $_[0] or 0;
	my $max = $_[1] or 0;
	my $decimals = $_[2] or 0;

	if($min eq "n")	{	$min = 0;	}
	if($max eq "n") {       $max = 0;       }
	if($decimals eq "n") {       $decimals = 0;       }

	my $range = $max-$min;
	my $num = rand($range);
	$num = sprintf("%.${decimals}f", $num+$min);
	return($num);
	} # end rand_between


sub rand_array
	{
	# Generates an array with a given number ($n) of random numbers between given values ($min, $max) and with a given number of decimals ($decimals)
	
	# Set error messages and accept input parameters
	my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
	my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$n, \$min, \$max, \$decimals)'\n\nwhere".
	"\t\$n is the number of numbers that should be returned. \$n cannot be 0.\n".
	"\t\$min is the lower limit of the number to be returned\n".
	"\t\$max is the upper limit of the number to be returned\n".
	"\t\$decimals is the number of desired decimals of the returned number\n\n";

	my @pars = @_ or die $usage;
	foreach my $el (@pars)  {       $el = text::trim($el);  }
	my $n = shift @pars or die $usage;      # NB! This will not work if the argument is the number 0 (because it will then be interpreted as false). In that case you need to use 'shift(@pars) or 0', but it may lead to proble$
	my $min = shift @pars or 0;
	my $max = shift @pars or 0;
	my $decimals = shift @pars or 0;

	# Reality check
	if($min >= $max)        {       die "The minimum value (${min}) is larger than the maximum value (${max}). Try again with the correct numbers!\n";  }

	# Initiate array to hold the results
	my @array=();

	# Generate numbers
	for(my $cc=1; $cc<=$n; $cc++)
		{
		my $random_number=stats::rand_between($min, $max, $decimals);
		push(@array, $random_number);
		}
	
	return(@array);
	} # end rand_array


sub mean_above
	{
	# Takes an array of numbers and checks if they have a mean above a given cutoff, and if the number of those numbers
	# that are above zero is larger than some given cutoff (e.g. the mean must be above 22.14 and the number of non-zero
	# samples must be at least 3). If the criteria are met the function returns "yes", else "no".

	my $subname="mean_above";
	my $usage="\nSyntax error for sub ${subname}. Correct usage: '\${packname}::\${subname}(\$arrayref, \$abund_limit, \$samp_limit)'\n\nwhere".
	"\t\$arrayref is a reference to an array holding the numbers whose mean you want to check is above some limit\n".
	"\t\$abund_limit is the minimum mean abundance for a feature to be accepted\n".
	"\t\$samp_limit is the minimum number of samples/libraries that need to have an abundance above 0 for a feature to be accepted\n\n";
	my @pars = @_ or die $usage;
	foreach my $el (@pars)  {       $el = text::trim($el);  }
	my $arref = shift @pars or die $usage;
	my $abund_limit = shift @pars or die $usage;
	my $samp_limit = shift @pars or die $usage;
	my @arr = @{$arref};

	my $total=0;
	my $above_0=0;
	for(my $cc=0; $cc<=$#arr; $cc++)
		{
		$total=$total+$arr[$cc];
		if($arr[$cc]>0)        {       $above_0++;     }
		}
	my $samp_num=scalar(@arr);
	my $mean=$total/$samp_num;
	if(($mean > $abund_limit) and ($above_0 >= $samp_limit))        {       return("yes");  }
	else    {       return("no");   }
	} # end mean_above


sub mean
	{
	# Calculates the arithmetic mean of an array of numbers
	
	# Set error messages and accept input parameters
	my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
	my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$numbers_array_ref)'\n\nwhere".
	"\t\$numbers_array_ref is a refernce to an array of numbers\n\n";
	my @pars = @_ or die $usage;
	foreach my $el (@pars)  {       $el = text::trim($el);  }
	my $arref = shift @pars or die $usage;
	my @arr = @{$arref};

	my $sum=0;
	my $n = scalar(@arr);
	foreach my $el (@arr)	{	$sum = $sum+$el;	}
	my $mean = $sum/$n;
	
	return($mean);
	} # end mean


sub num_above_cutoff
	{
	# Takes an array of numerical values and a numerical cutoff value, and returns the number of values that are
	# larger than or equal to that cutoff, as well as the percentage of values larger than or equal to the cutoff.

	# Set error messages and accept input parameters
	my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
	my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$numbers_array_ref, \$cutoff)'\n\nwhere".
	"\t\$numbers_array_ref is a refernce to an array of numbers\n".
	"\t\$cutoff is the value, e.g. 0, that a number must be above to count as being above that value\n\n";
	my @pars = @_ or die $usage;
	foreach my $el (@pars)	{	$el=text::trim($el);	}
	my $numbers_array_ref = shift(@pars) or die $usage;
	my $cutoff = shift(@pars) or 0;

	my @arr = @{$numbers_array_ref};
	my $num_above_cutoff=0;

	for(my $cc=0; $cc<=$#arr; $cc++)
		{
		if($arr[$cc]>=$cutoff)	{	$num_above_cutoff++;	}		
		}

	my $perc_above_cutoff = 100*($num_above_cutoff/scalar(@arr));
	my @stats=($num_above_cutoff, $perc_above_cutoff);

	return(@stats);
	} # end num_above_cutoff


sub median
	{
	# Computes the median of a set of numbers
	
	# Set error messages and accept input parameters
	my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
	my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$arrayref)'\n\nwhere".
	"\t\$arrayref is a reference to an array of numbers whose median should be computed\n\n";
	my @pars = @_ or die $usage;
	foreach my $el (@pars)  {       $el = text::trim($el);  }
	my $arref = shift @pars or die $usage;      # NB! This will not work if the argument is the number 0 (because it will then be interpreted as false). In that case you need to use 'shift(@pars) or 0', but it may lead to proble$

	my @arr=@{$arref};
	my $num_obs=scalar(@arr);
	@arr = sort { $a <=> $b } @arr;
	my $median_index = ($num_obs/2)-1;	# -1 is because perl starts counting aray indices at 0 rather than 1
	my $median=0;

	# If the number of observations is even
	if(int($median_index) == $median_index)	{	$median = ($arr[$median_index]+$arr[1+$median_index])/2;	}

	# If the number of observations is odd	
	elsif(int($median_index) != $median_index)
		{
		$median_index = int($median_index)+1;
		$median = $arr[$median_index];
		}

	return($median);
	} # end median


sub min
	{
	# Returns the minimum value of a set of numbers
	
	# Set error messages and accept input parameters
	my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
	my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$arrayref)'\n\nwhere".
	"\t\$arrayref is a reference to an array of numbers whose minimum value should be returned\n\n";
	my @pars = @_ or die $usage;
	foreach my $el (@pars)  {       $el = text::trim($el);  }
	my $arref = shift @pars or die $usage;      # NB! This will not work if the argument is the number 0 (because it will then be interpreted as false). In that case you need to use 'shift(@pars) or 0', but it may lead to proble$
	
	my @arr = @{$arref};

	# Sort the array numerically
	@arr = sort { $a <=> $b } @arr;

	# Get the minimum value
	my $min = shift(@arr);

	return($min);
	} # end min


sub max
	{
	# Returns the maximum value of a set of numbers
	
	# Set error messages and accept input parameters
	my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
	my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$arrayref)'\n\nwhere".
	"\t\$arrayref is a reference to an array of numbers whose maximum value should be returned\n\n";
	my @pars = @_ or die $usage;
	foreach my $el (@pars)  {       $el = text::trim($el);  }
	my $arref = shift @pars or die $usage;      # NB! This will not work if the argument is the number 0 (because it will then be interpreted as false). In that case you need to use 'shift(@pars) or 0', but it may lead to proble$

	my @arr = @{$arref};

	# Sort the array numerically
	@arr = sort { $b <=> $a } @arr;

	# Get the maximum value
	my $max = shift(@arr);

	return($max);
	} # end max


sub variance
	{
	# Returns the variance of a set of numbers, specified as an array reference. Also specify whether to compute
	# the sample variance (set $samp_or_pop = “samp”) or the population variance ( $samp_or_pop = “pop”).

	# Set error messages and accept input parameters
	my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
	my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$arrayref, \$samp_or_pop)'\n\nwhere".
	"\t\$arrayref is a reference to an array of numbers whose variance should be returned\n".
	"\t\$samp_or_pop is an indicator of whether the sample or population variance should be computed\n\n";
	my @pars = @_ or die $usage;
	foreach my $el (@pars)  {       $el = text::trim($el);  }
	my $arref = shift @pars or die $usage;
	my $samp_or_pop = shift @pars or die $usage;
	
	my @arr = @{$arref};
	my $mean=stats::mean($arref);
	my $sum_of_squares=0;

	for(my $cc=0; $cc<=$#arr; $cc++)
		{
		my $diff = $arr[$cc]-$mean;
		my $square = $diff ** 2;
		$sum_of_squares += $square;
		}

	my $n=0;
	if($samp_or_pop eq "samp")	{	$n = (scalar(@arr))-1;	}
	elsif($samp_or_pop eq "pop")	{	$n = scalar(@arr);	}
	else
		{
		die "Parameter $samp_or_pop for subroutine $subname (called by script '${calling_script}', line ${calling_line}) ".
		"needs to be either 'samp' or 'pop'. Check that it is!\n";
		}
	my $variance=0;
	unless($n==0)	{	$variance = $sum_of_squares/$n;	}
	return($variance);
	} # end variance


sub stdev
	{
	# Returns the standard deviation of a set of numbers
	
	# Set error messages and accept input parameters
	my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
	my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$arrayref, \$samp_or_pop)'\n\nwhere".
	"\t\$arrayref is a reference to an array of numbers whose standard deviation should be returned\n".
	"\t\$samp_or_pop is an indicator of whether the sample or population standard deviation should be computed\n\n";
	my @pars = @_ or die $usage;
	foreach my $el (@pars)  {	$el = text::trim($el);  }
	my $arref = shift @pars or die $usage; 
	my $samp_or_pop = shift @pars or die $usage;
	unless(($samp_or_pop eq "samp") or ($samp_or_pop eq "pop"))
		{
		die "Parameter $samp_or_pop for subroutine $subname (called by script '${calling_script}', line ${calling_line}) ".
		"needs to be either 'samp' or 'pop'. Check that it is!\n";
		}

	my @arr = @{$arref};
	my $variance=stats::variance($arref, $samp_or_pop);
	my $stdev = sqrt($variance);

	return($stdev);
	} # end stdev


sub range_string
	{
	# Returns a string representing the range of a set of numbers in the format "start_stop"
	
	# Set error messages and accept input parameters
	my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
	my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$arrayref)'\n\nwhere".
	"\t\$arrayref is a reference to an array of numbers whose range should be returned\n\n";
	my @pars = @_ or die $usage;
	foreach my $el (@pars)  {       $el = text::trim($el);  }
	my $arref = shift @pars or die $usage;

	my @arr = @{$arref};
	my $min=stats::min($arref);
	my $max=stats::max($arref);
	my $range = $min."_".$max;

	return($range);
	} # end range_string


sub percentile
	{
	# Computes the k:th percentile of a set of numbers
	
	# Set error messages and accept input parameters
	my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
	my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$arrayref, \$k)'\n\nwhere".
	"\t\$arrayref is a reference to an array holding the set of numers you want to compute a percentile for\n".
	"\t\$k is the percentile to be computed, e.g. 75\n\n";

	my @pars = @_ or die $usage;
	foreach my $el (@pars)  {       $el = text::trim($el);  }
	my $arref = shift @pars or die $usage;      # NB! This will not work if the argument is the number 0 (because it will then be interpreted as false). In that case you need to use 'shift(@pars) or 0', but it may lead to proble$
	my $k = shift @pars or die $usage;
	my $fraction = $k/100;

	my @arr = @{$arref};
	@arr = sort { $a <=> $b} @arr;
	
	# How many values are there?
	my $n = scalar(@arr);

	# Which value should I pick?
	my $index = $fraction * $n;
	my $whole = int($index);
	my $decimal = $index-$whole;

	# If index is not a whole number, round it
	if($index != $whole)
		{
		if($decimal >= 0.5)	{	$index = $index+1;	}
		}

	# Get percentile
	my $percentile = $arr[$index];

	return($percentile);	
	} # end percentile


sub ofc
	{
	# Compares two sets of numbers by computing an offset fold change bewteen them
	
	# Set error messages and accept input parameters
	my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
	my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$array_reference_1, \$array_reference_2, \$offset)'\n\nwhere".
	"\t\$array_reference_1 is a reference to an array holding the first set of numbers that should be compared\n".
	"\t\$array_reference_2 is a reference to an array holding the second set of numbers that should be compared\n".
	"\t\$offset is the offset to be used\n\n";
	my @pars = @_ or die $usage;
	foreach my $el (@pars)  {       $el = text::trim($el);  }
	my $arref1 = shift @pars or die $usage;      # NB! This will not work if the argument is the number 0 (because it will then be interpreted as false). In that case you need to use 'shift(@pars) or 0', but it may lead to proble$
	my $arref2 = shift @pars or die $usage;
	my $offset = shift @pars or die $usage;

	my $min_1=stats::min($arref1);
	my $max_1=stats::max($arref1);
	my $min_2=stats::min($arref2);
	my $max_2=stats::max($arref2);

	my $ofc="Faster_Tinne";

	# If range 1 is higher than range 2
	if($min_1>$max_2)
		{	
		my $ratio = ($offset+$min_1)/($offset+$max_2);		
		$ofc = (log($ratio))/(log(2));
		}

	# If range 2 is higher than range 1
	elsif($min_2>$max_1)
		{
		my $ratio = ($offset+$min_2)/($offset+$max_1);
		my $ofc2 = (log($ratio))/(log(2));
		$ofc = -$ofc2;
		}
	else	{	$ofc = 0;	}
	return($ofc);
	} # end ofc


sub list_overlap
	{
	# Compares two lists and returns the elements they have in common
	
	# Set error messages and accept input parameters
	my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
	my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$array_reference_1, \$array_reference_2)'\n\nwhere".
	"\t\$array_reference_1 is a reference to an array holding the first set of numbers that should be compared\n".
	"\t\$array_reference_2 is a reference to an array holding the second set of numbers that should be compared\n";
	my @pars = @_ or die $usage;
	foreach my $el (@pars)  {       $el = text::trim($el);  }
	my $arref1 = shift @pars or die $usage;      # NB! This will not work if the argument is the number 0 (because it will then be interpreted as false). In that case you need to use 'shift(@pars) or 0', but it may lead to proble$
	my $arref2 = shift @pars or die $usage;

	my @arr1 = @{$arref1};
	my @arr2 = @{$arref2};
	
	my @common=();	

	for(my $c=0; $c<=$#arr1; $c++)
		{
		for(my $d=0; $d<=$#arr2; $d++)
			{
			if($arr1[$c] eq $arr2[$d])	{	push(@common, $arr1[$c]);	}
			}
		}

	return(@common);
	} # end list_overlap


sub test_diff_expression_range
	{
	# Tests differential expression of genes between groups of samples in an RNAseq experiment, based on non-overlap of the ranges of expression
	# for individual genes between groups of samples. Takes as input a normalised expression matrix in csv format and a csv file listing each sample
	# to be included in the ananlysis and the groups to which these samples belong. It also needs to list the group-group comparisons to be made. Example
	# of required format:
	#
	#        sample1,groupA
	#        sample2,groupA
	#        sample3,groupB
	#        sample4,groupB
	#        sample5,groupC
	#        sample6,groupC
	#        groupA_vs_groupB
	#        groupA_vs_groupC
	#        groupB_vs_groupC
	#
	# Outputs a new expression matrix with columns for offset fold change (OFC) and significance of OFC for each group-group comparison added.
	
	# Set error messages and accept input parameters
	my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
	my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$expression_matrix.csv, \$groups.csv, \$offset, \$header_y_n, \$outname)'\n\nwhere".
	"\t\$expression_matrix.csv is a gene expression matrix in csv format. There needs to be a header line that includes sample names\n".
	"\t\$groups.csv is a csv file listing each sample to be included in the ananlysis, the groups to which thise samples belong, and which group comparisons should be made, e.g\n".
	"\tsample_1,group_A\n".
	"\tsample_2,group_A\n".
	"\tsample_3,group_B\n".
	"\tsample_4,group_B\n".
	"\tsample_5,group_C\n".
	"\tsample_6,group_C\n".
	"\tgroupA_vs_groupB\n".
	"\tgroupA_vs_groupC\n".
	"\tgroupB_vs_groupC\n\n".
	"\t\$offset is the offset to be used in the computation of foldchange. If this is to be 0, set it to 'n' (perl can't handle zeros as command line arguments)\n".
	"\t\$header_y_n is an indicator of whether the input expression matrix has headers or not. If not, sample names and order will be taken from the groups.csv file\n".
	"\t\$outname is the desired name of the outfile, which is a new expression matrix with two added columns (OFC and significance) for each group-group comparison made\n\n";
	my @pars = @_ or die $usage;
	foreach my $el (@pars)  {       $el = text::trim($el);  }
	my $infile = shift @pars or die $usage;      # NB! This will not work if the argument is the number 0 (because it will then be inte	rpreted as false). In that case you need to use 'shift(@pars) or 0', but it may lead to proble$
	my $group_file = shift @pars or die $usage;
	my $offset = shift @pars or die $usage;
	my $header_y_n = shift @pars or die $usage;
	my $outname = shift @pars or die $usage;

	if($offset eq "n")	{	$offset=0;	}
	my $splice_here=0;

	# Read the sample names file
	my ($group_matrix_ref, $comparisons_ref)=fileTools::read_diff_expr_sample_list($group_file);
	my @group_matrix = @{$group_matrix_ref};
	my @comparisons_matrix = @{$comparisons_ref};	

	# Read the sample names file again to extract header information
	my @header_matrix=fileTools::read_table($group_file, "csv");
	my @header=();
	push(@header, "Feature_ID");
	for(my $zz=0; $zz<=$#header_matrix; $zz++)
		{
		if($header_matrix[$zz][0] !~ /_vs_/)	{	push(@header, $header_matrix[$zz][0]);	}
		}
	
	# Read the expression matrix
	my @matrix=fileTools::read_table($infile, "csv");
	if($header_y_n eq "y")
		{
		my $header_ref = shift(@matrix);	
		@header = @{$header_ref};
		}

	my @new_header = @header;
	my @new_matrix=();

	# Search the @header for each sample name and get the indices (column numbers) for each group of samples
	my @group_indices=();

	# For each group...
	for(my $cc=0; $cc<=$#group_matrix; $cc++)
		{
		my @arr = @{$group_matrix[$cc]};
		$group_indices[$cc][0]=$group_matrix[$cc][0];

		# For each sample in group...
		for(my $dd=1; $dd<=$#arr; $dd++)
			{
			my $samp_name=$arr[$dd];
			# Find the column number of the current sample
			my ($samp_ind) = grep { $header[$_] eq "$samp_name" } 0..$#header;
			if($samp_ind > $splice_here)	{	$splice_here = $samp_ind;	}			
			$group_indices[$cc][$dd]=$samp_ind;
			}
		}

	$splice_here++;
	my $splice_start = $splice_here;

	# Loop over genes in matrix
	for(my $ee=0; $ee<=$#matrix; $ee++)
		{
		my @arr = @{$matrix[$ee]};

		# For each comparison to be made...
	        for(my $zz=0; $zz<=$#comparisons_matrix; $zz++)
			{
			# Get group names
			my $first_group_name = $comparisons_matrix[$zz][0];
			my $second_group_name = $comparisons_matrix[$zz][1];
			my @first_group_indices=();
			my @second_group_indices=();
			my @first_group_data=();
			my @second_group_data=();
			my $comp_sig="NO";

			# Get data for the two groups
			for(my $ff=0; $ff<=$#group_indices; $ff++)
				{
				my @group = @{$group_indices[$ff]};
				if($group[0] eq $first_group_name)
					{
					@first_group_indices=@group[1..$#group];	
					@first_group_data=@arr[@first_group_indices];					
					}
				elsif($group[0] eq $second_group_name)
					{
					@second_group_indices=@group[1..$#group];
					@second_group_data=@arr[@second_group_indices];
					}
				}
			
			# Compute offset foldchange
			my $ofc=stats::ofc(\@first_group_data, \@second_group_data, $offset);
			if($ofc != 0)	{	$comp_sig = "YES";	}

			# Splice the results into the current gene array (and, if this is the first gene, also column descriptions into @header)
			splice(@arr, $splice_here, 0, $comp_sig, $ofc);
			if($ee==0)	{	splice(@new_header, $splice_here, 0, "${first_group_name}_${second_group_name}_sign", "${first_group_name}_${second_group_name}_OFC");	}
			$splice_here += 2;

			} # End of 'for each comparison to be made'	
	
		$splice_here = $splice_start;		
		push(@new_matrix, [@arr]);

        } # Loop over genes ends

	# Add @new_header to @new_matrix
	unshift(@new_matrix, \@new_header);

	# Print the new matrix to an outfile
	fileTools::write_table(\@new_matrix, "csv", $outname, "lin");
	} # end test_diff_expression_range


sub geometric_mean
	{
	# Returns the geometric mean of an array of numbers
	
	# Set error messages and accept input parameters
	my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
	my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$array_ref)'\n\nwhere".
	"\t\$array_ref is a refernce to the array holding the numbers for which a geometric mean should be computed\n\n";
	my @pars = @_ or die $usage;
	foreach my $el (@pars)  {       $el = text::trim($el);  }
	my $arref = shift @pars or die $usage;      # NB! This will not work if the argument is the number 0 (because it will then be interpreted as false). In that case you need to use 'shift(@pars) or 0', but it may lead to proble$

	my @arr = @{$arref};
	my $n = scalar(@arr);
	my $product=1;
	for(my $cc=0; $cc<=$#arr; $cc++)	{	$product = $product*$arr[$cc];	}
	my $geomean = $product**(1/$n);
	return($geomean);
	} # end geometric_mean


sub three_point_expr_pattern
	{
	# Takes an analysed expression matrix and and adds a column for three-point expression pattern class. Requires the matrix to have a header with, for each
	# differential expression comparison (bewteen treatment groups/conditions) made, two columns; one for statistical significance and one for offset fold change.
	# The labels of these columns must contain the text "_sign" and "_OFC" respectively (that is how the function knows which values to use when looking for patterns).

	# Set error messages and accept input parameters
	my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
	my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$infile.csv, \$ofc_cutoff, \$outname)'\n\nwhere".
	"\t\$infile.csv is the matrix to be analysed\n". 
	"\t\$ofc_cutoff is the offset foldchange value that differences must be greater than to be considered of interest\n". 
	"\t\$outname is the preferred name of the outfile\n\n";
	my @pars = @_ or die $usage;
	foreach my $el (@pars)  {       $el = text::trim($el);  }
	my $infile = shift @pars or die $usage;      # NB! This will not work if the argument is the number 0 (because it will then be interpreted as false). In that case you need to use 'shift(@pars) or 0', but it may lead to proble$
	my $ofc_cutoff = shift @pars or die $usage;
	my $outname = shift @pars or die $usage;
	
	# Read matrix
	my @matrix=fileTools::read_table($infile, "csv");
	
	my $headerref = shift(@matrix);
	my @header = @{$headerref};
	my $last_ind = scalar(@header);
	push(@header, "pattern_class", "pattern_class_name");

	# Find significance and ofc columns
	my (@sig_inds) = grep { $header[$_] =~ /_sign/ } 0..$#header;
	my (@ofc_inds) = grep { $header[$_] =~ /_OFC/ } 0..$#header;

	#my $last_ind="Faster_Tinne";

	# Loop over matrix
	for(my $cc=0; $cc<=$#matrix; $cc++)
		{
		my @arr = @{$matrix[$cc]};
		#if($cc==0)	{	$last_ind = scalar(@arr);	}
		my $class_name=1;
		my $class="Faster_Tinne";
		my $sig_first=$matrix[$cc][$sig_inds[0]];
		my $sig_second=$matrix[$cc][$sig_inds[1]];
		my $ofc_first=$matrix[$cc][$ofc_inds[0]];
		my $ofc_second=$matrix[$cc][$ofc_inds[1]];
		
		# Evaluate conditions and assign class

		# If first is the same
		if($sig_first eq "NO")
			{
			# and second is the same
			if($sig_second eq "NO")	{	$class=1; $class_name="same-same";	}
			
			# and second is not the same
			elsif($sig_second eq "YES")
				{
				# If second is positive
				if($ofc_second >= $ofc_cutoff)	{	$class=2; $class_name="same-up";}
				# If second is negative
				elsif($ofc_second <= -$ofc_cutoff)	{	$class=3; $class_name="same-down";	}
				# If second is too small
				else	{	$class=1; $class_name="same-same";	}
				}
			}

		# If first is not the same
		elsif($sig_first eq "YES")
			{
			# If first is positive
			if($ofc_first >= $ofc_cutoff)
				{
				# If second is the same
				if($sig_second eq "NO")	{	$class=4; $class_name="up-same";	}

				# If second is not the same				
				elsif($sig_second eq "YES")
					{
					# If second is positive
					if($ofc_second >= $ofc_cutoff)	{	$class=5; $class_name="up-up";	}
					# If second is negative
					elsif($ofc_second <= -$ofc_cutoff)	{	$class=6; $class_name="up-down";	}
					# If second is too small
					else	{	$class=4; $class_name="up-same";	}
					}
				}

			# If first is negative
			elsif($ofc_first <= -$ofc_cutoff)
				{
				# If second is the same
				if($sig_second eq "NO") {       $class=7; $class_name="down-same";        }

				# If second is not the same
				elsif($sig_second eq "YES")
					{
					# If second is positive
					if($ofc_second >= $ofc_cutoff)  {       $class=8; $class_name="down-up";  }
					# If second is negative
					elsif($ofc_second <= -$ofc_cutoff)      {       $class=9; $class_name="down-down";        }
					# If second is too small
					else    {       $class=7; $class_name="down-same";        }
					}
				}

			# If first is too small
			else
				{
				# If second is the same
				if($sig_second eq "NO") {       $class=1; $class_name="same-same";        }

				# If second is not the same
				elsif($sig_second eq "YES")
					{
					# If second is positive
					if($ofc_second >= $ofc_cutoff)  {       $class=2; $class_name="same-up";  }
					# If second is negative
					elsif($ofc_second <= -$ofc_cutoff)      {       $class=3; $class_name="same-down";        }
					# If second is too small
					else    {       $class=1; $class_name="same-same";        }
					}
				}
			}

		$matrix[$cc][$last_ind]=$class;
		$matrix[$cc][1+$last_ind]=$class_name;
		}	
	
	unshift(@matrix, \@header);
	fileTools::write_table(\@matrix, "csv", $outname, "lin");
    } # end three_point_expr_pattern


sub assign_ranks
	{
	# This function takes (a reference to) an array of values, sorts them in ascending order and assigns ranks to them. If there are tied values in
	# the input array, the mean rank of the values in a group of tied values can be assigned to those values. Set argument $adjust_ties_y_n to 'y'
	# if you want this, or to 'n' if you want to use the original , unadjusted ranks. The function returns an array of the ranks in the original order
	# of elements in the input array. If your array has numeric values, set argument $num_or_char to "num", otherwise to "char".

	# Set error messages
	my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
	my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$array_ref, \$adjust_ties_y_n, \$num_or_char)'\n\nwhere".
	"\t\$array_ref is a reference to an array holding the numbers or strings to be ranked\n".
	"\t\$adjust_ties_y_n  Set this to yes if you want values in a group of ties (=equal) values to be assigned the average rank within that group\n".
	"\t\$num_or_char  Set this to 'num' if the values to be ranked are numeric, or 'char' if they are text\n\n";
	
	# Accept input parameters
	my @pars = @_ or die $usage;
	foreach my $el (@pars)  {       $el = text::trim($el);  }
	my $arref = shift @pars or die $usage;
	my $adj_ties = shift @pars or die $usage;
	my $num_or_char = shift @pars or die $usage;
	my @arr = @{$arref};

	# Make an array with the original order of values in @arr, and combine @orig_order and @arr into a matrix
	my @orig_order = (1..scalar(@arr));	
	my @matrix=();
	
	@matrix=matrixTools::add_cols_to_matrix(\@matrix, \@orig_order, \@arr);

	# Sort matrix ascending on @arr
	if($num_or_char eq "num")	{	@matrix = sort { $a->[1] <=> $b->[1] } @matrix;	}
	elsif($num_or_char eq "char")	{	@matrix = sort { $a->[1] cmp $b->[1] } @matrix;	}
	
	# Assign preliminary ranks to values in @arr (this will be the third column of @matrix)
	my @prel_ranks = (1..scalar(@matrix));
	@matrix=matrixTools::add_cols_to_matrix(\@matrix, \@prel_ranks);

	my @outranks=();
	if($adj_ties eq "n")
		{
		# Sort matrix back into original order
		@matrix = sort { $a->[0] <=> $b->[0] } @matrix;
		@outranks = matrixTools::get_matrix_columns_as_rows(\@matrix, 3);
		}
	
	elsif($adj_ties eq "y")
		{
		# Loop over matrix to assign correct ranks to tied values (=equal values)
		my @tied_groups=();		# This is a matrix that will hold arrays. Each array is the ranks of a group of tied values.
		my @temp_group=();
		for(my $cc=0; $cc<=$#matrix; $cc++)
			{
			# If either the value before or the value after is equal, call this a tied value,
			# and add it's index to @temp_group
			
			# If this is the first line of @matrix, compare only to the next value (the previous one isn't defined)
			if($cc==0)
				{
				# If this is a tied value...
				if($num_or_char eq "num")	{	if($matrix[$cc][1] == $matrix[${cc}+1][1])	{	push(@temp_group, $cc);	}	}
				elsif($num_or_char eq "char")	{	if($matrix[$cc][1] eq $matrix[${cc}+1][1])	{	push(@temp_group, $cc);	}	}
				} # "If first line of matrix" ends
			
			# If this is the last line of @matrix, compare only to the previous value (the next one isn't defined)
			elsif($cc==$#matrix)
				{
				if($num_or_char eq "num")
					{
					if($matrix[$cc][1] == $matrix[${cc}-1][1])
						{
						push(@temp_group, $cc);
						push(@tied_groups, [@temp_group]);
						@temp_group=();
						}
					}
					
				elsif($num_or_char eq "char")
					{
					if($matrix[$cc][1] eq $matrix[${cc}-1][1])
						{
						push(@temp_group, $cc);
						push(@tied_groups, [@temp_group]);
						@temp_group=();
						}
					}
				} # "If last line of @matrix" ends here
				
			# If it is any other line
			else
				{
				if($num_or_char eq "num")
					{
					if(($matrix[$cc][1] == $matrix[${cc}-1][1]) or ($matrix[$cc][1] == $matrix[${cc}+1][1]))
						{
						push(@temp_group, $cc);
						
						# If the value after this one isn't the same, end tied group
						if($matrix[$cc][1] != $matrix[${cc}+1][1])
							{
							push(@tied_groups, [@temp_group]);
							@temp_group=();
							}
						} # "If value is tied" ends here	
					}
					
				elsif($num_or_char eq "char")
					{
					if(($matrix[$cc][1] eq $matrix[${cc}-1][1]) or ($matrix[$cc][1] eq $matrix[${cc}+1][1]))
						{
						push(@temp_group, $cc);
						
						# If the value after this one isn't the same, end tied group
						if($matrix[$cc][1] ne $matrix[${cc}+1][1])
							{
							push(@tied_groups, [@temp_group]);
							@temp_group=();
							}
						} # "If value is tied" ends here	
					}					
				} # "If any other line" ends here
			} # Loop over matrix rows ends here
		
		# We now have a @tied_groups matrix whose rows are lists of the @matrix rows that have tied values.
		
		# Loop over rows in @tied_groups matrix
		for(my $dd=0; $dd<=$#tied_groups; $dd++)
			{
			my @group_indices = @{$tied_groups[$dd]};
			my @group=();
		
			# Loop over indices in the current group of tied ranks and compute the mean rank
			for(my $ee=0; $ee<=$#group_indices; $ee++)
				{
				my $matrix_row = $group_indices[$ee];
				my $tied_rank = $matrix[$matrix_row][2]; 
				push(@group, $tied_rank);
				}
			my $mean_rank=stats::mean(\@group);
			
			# Loop over indices in the current group of tied ranks again, and assign the mean rank to all values in @matrix that belong to the group
			for(my $ff=0; $ff<=$#group_indices; $ff++)
				{
				my $matrix_row2 = $group_indices[$ff];
				$matrix[$matrix_row2][3] = $mean_rank;
				}
			}

		# Loop over @matrix again and set adjusted rank for all values that are not tied to the preliminary rank for that value
		for(my $gg=0; $gg<=$#matrix; $gg++)
			{
			if(!defined($matrix[$gg][3]))	{	$matrix[$gg][3] = $matrix[$gg][2];	}
			}
			
		# Sort matrix back to original order
		if($num_or_char eq "num")	{	@matrix = sort { $a->[0] <=> $b->[0] } @matrix;	}
		elsif($num_or_char eq "char")	{	@matrix = sort { $a->[0] cmp $b->[0] } @matrix;	}
		
		# Get the adjusted ranks column from matrix
		my @adj_ranks=matrixTools::get_matrix_columns_as_rows(\@matrix, 4);
		@outranks = @adj_ranks;
		}
	
	return(@outranks);
	} # end assign_ranks


sub covariance
	{
	# Computes the covariance between two sets of numbers (specified as array references). The two sets must have the same number of numbers.
	# Specify whether the sample or population covariance should be computed (set $samp_or_pop to ‘samp’ or ‘pop’ respectively).

	# Set error messages
	my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
	my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$arrayref1, \$arrayref2, \$samp_or_pop)'\n\nwhere".
	"\t\$arrayref1, \$arrayref2 are references to the two arrays (lists) holding the numbers to be used in the computation\n".
	"\t\$samp_or_pop is an indicator of whether the sample (set to 'samp') or population (set to 'pop') covariance should be computed\n\n";
	
	# Accept input parameters
	my @pars = @_ or die $usage;
	foreach my $el (@pars)  {       $el = text::trim($el);  }
	my $arref1 = shift @pars or die $usage;
	my $arref2 = shift @pars or die $usage;
	my $samp_or_pop = shift @pars or die $usage;

	my @arr1 = @{$arref1};
	my @arr2 = @{$arref2};
	
	# Get the means of the two sets of numbers
	my $mean1=stats::mean($arref1);
	my $mean2=stats::mean($arref2);
	
	# Compute the degrees of freedom
	my $df=0;
	if($samp_or_pop eq "samp")	{	$df = (scalar(@arr1))-1;	}
	elsif($samp_or_pop eq "pop")	{	$df = scalar(@arr1);	}
	else
		{
		die "Parameter \$samp_or_pop for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}) ".
		"needs to be either 'samp' or 'pop'. Please check that it is.\n"
		}

	# Loop over numbers in the two sets and compute products of deviations
	my $sum_of_products=0;
	
	for(my $cc=0; $cc<=$#arr1; $cc++)
		{
		my $dev1 = $arr1[$cc] - $mean1;
		my $dev2 = $arr2[$cc] - $mean2;
		my $product = $dev1 * $dev2;
		$sum_of_products = $sum_of_products + $product;
		}
	
	# Divide the sum of products by the degrees of freedom
	my $covariance = $sum_of_products/$df;
	
	return($covariance);
	} # end covariance
	

sub corr_coeff
	{
	# Computes the correlation coefficient between two sets of numbers (specified as array references). The two sets must have the same number of numbers.
	# Specify whether the sample or population correlation coefficient should be computed (set $samp_or_pop to ‘samp’ or ‘pop’ respectively).
	# Currently, the Pearson Product Moment Correlation Coefficient is the one returned, but more options (e.g. Spearman) should be added in the future.

	# Set error messages
	my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
	my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$arrayref1, \$arrayref2, \$samp_or_pop)'\n\nwhere".
	"\t\$arrayref1, \$arrayref2 are references to the two arrays (lists) holding the numbers to be used in the computation\n".
	"\t\$samp_or_pop is an indicator of whether the sample (set to 'samp') or population (set to 'pop') correlation coefficient should be computed\n\n";
	
	# Accept input parameters
	my @pars = @_ or die $usage;
	foreach my $el (@pars)  {       $el = text::trim($el);  }
	my $arref1 = shift @pars or die $usage;
	my $arref2 = shift @pars or die $usage;
	my $samp_or_pop = shift @pars or die $usage;
	#my $type = shift @pars or die $usage;
	my $type = "pmcc";

	my @arr1 = @{$arref1};
	my @arr2 = @{$arref2};
	my $corr_coeff = 0;
	
	if($type eq "pmcc")
		{
		my $stdev1=stats::stdev($arref1, $samp_or_pop);
		my $stdev2=stats::stdev($arref2, $samp_or_pop);
		my $denominator = $stdev1 * $stdev2;
		my $covariance=stats::covariance($arref1, $arref2, $samp_or_pop);
		
		print "stdev1 = $stdev1\n";
		print "stdev2 = $stdev2\n";
		print "stdev product = $denominator\n";
		print "covariance = $covariance\n";
		
		$corr_coeff = $covariance/$denominator;
		}
	
	return($corr_coeff);
	} # end corr_coeff
	
sub pearson_corr
	{
	# Computes the Pearson Product Moment Correlation Coefficient between two sets of numbers (specified as array references). The two sets must have
	# the same number of numbers. Specify whether the sample or population correlation coefficient should be computed (set $samp_or_pop to ‘samp’
	# or ‘pop’ respectively).

	# Set error messages
	my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
	my $usage1="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}).\n";
	my $usage2 = "Correct usage: '${subname}(\$arrayref1, \$arrayref2, \$samp_or_pop)'\n\nwhere".
	"\t\$arrayref1, \$arrayref2 are references to the two arrays (lists) holding the numbers to be used in the computation\n".
	"\t\$samp_or_pop is an indicator of whether the sample (set to 'samp') or population (set to 'pop') correlation coefficient should be computed\n\n";
	
	# Accept input parameters
	my @pars = @_ or die $usage1, $usage2;
	foreach my $el (@pars)  {       $el = text::trim($el);  }
	my $arref1 = shift @pars or die $usage1, $usage2;
	my $arref2 = shift @pars or die $usage1, $usage2;
	my $samp_or_pop = shift @pars or die $usage1, $usage2;
	my @arr1 = @{$arref1};
	my @arr2 = @{$arref2};
	
	# Assert that input parameters are valid
	if(scalar(@arr1) != scalar(@arr2))	{	die $usage1, "The two input arrays must be of equal length\n";	}
	unless(listTools::value_in_array($samp_or_pop, ["samp", "pop"], "char"))	{	die $usage1, "Parameter \$samp_or_pop must be either 'samp' or 'pop'\n";	}

	my $corr_coeff = 0;
	
	my $stdev1=stats::stdev($arref1, $samp_or_pop);
	my $stdev2=stats::stdev($arref2, $samp_or_pop);
	my $denominator = $stdev1 * $stdev2;
	my $covariance=stats::covariance($arref1, $arref2, $samp_or_pop);
	
	$corr_coeff = $covariance/$denominator;
	
	print "\n\t$corr_coeff\n";
			
	return($corr_coeff);
	} # end pearson_corr


return(1);

# end functions