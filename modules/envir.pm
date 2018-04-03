package envir;

# timestamp()
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


sub timestamp
	{
	# Create a timestamp
	my ($sec,$min,$hour,$mday,$month,$year,$wday,$yday,$isdst) = localtime(time);
	$year=$year+1900;
	$month=$month+1;
	if($month<10)   {       $month = "0"."$month";  }
	if($mday<10)   {       $mday = "0"."$mday";  }
	my $timestamp = "$year"."-"."$month"."-"."$mday"."_"."$hour"."-"."$min"."-"."$sec";

	return($timestamp);
	} # end timestamp
	
sub check_arg
	{
	# Validates user input to scripts and functions. Returns false if everything is ok, and an error message otherwise.
	# (The decision to terminate processing lies in the script that calls this subroutine, rather than in the subroutine
	# itself).
	
	# Usage: subname($argument, [allowed_type1, allowed_type2, ...], [allowed_value1, allowed_value2, ...]);
	
	# Input parameters:
	
	# $argument	=	the argument that should be validated. If argument is 0, set to 'zero'
	# allowed_type		=	the argument's type. Options: "arrayref", "matrixref", "string", "number"
	# allowed_value		=	a value that the argument can take. If not relevant, leave blank
	
	
	# Set usage message
	my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
	my $sub_and_caller = "'${subname}' (called by script '${calling_script}', line ${calling_line})";
	my $usage="\nUsage error for subroutine $sub_and_caller\n".
	"Correct usage: '${subname}(\$argument, [allowed_type1, allowed_type2,...], [allowed_value1, allowed_value2,...])'\n\nwhere".
	
	"\targument\t=\tthe argument that should be validated. If argument is 0, set to 'zero'\n".
	"\tallowed_type\t=\tthe argument's type. Options: 'char', 'num', 'arrayref', 'matrixref', 'hashref'\n".
	"\tallowed_value\t=\ta value that the argument can take. If not relevant, leave blank\n\n)";

	# Read input parameters
	my @pars = @_ or die $usage;
	foreach my $el (@pars)  {       $el = text::trim($el);  }
	
	my $argument = shift @pars or die $usage;
	my $type_options_ref = shift @pars or die $usage;
	my $value_options_ref = shift @pars or 0;
	
	# Process parameters
	my $message = "";
	my @type_options = @{$type_options_ref};
	my @value_options;
	if($value_options_ref)	{	@value_options = @{$value_options_ref};	}		# If the user has supplied allowed values, put them in an array
	
	my $dec=`perl -MPOSIX=locale_h -e "print localeconv()->{decimal_point}`;	# Get the decimal separator for the current system
	
	# Check that the argument's type is one those allowed
	my $argtype=misc::type($argument);
	if($argument eq "zero")	{	$argtype = "num";	}
	if(not listTools::value_in_array($argtype, $type_options_ref, "char"))	{	$message = "The argument has type $argtype, but the allowed types are (@type_options)";	}
	
	# If type is ok, check allowed values
	else
		{
		# Check that the argument's value is one of those allowed (if specified)
		if($value_options_ref)
			{
			# If the argument is a string
			if($argtype eq "char")
				{
				unless(listTools::value_in_array($argument, $value_options_ref, "char"))	{	$message = "The argument has value \"$argument\", but the allowed values are \"@value_options\"";	}
				}
				
			# Else, if the argument is numeric
			elsif($argtype eq "num")
				{
				if($argument eq "zero")	{	$argument = 0;	}
				
				# If the allowed values are equalities/inequalities
				my $equality;
				foreach my $element (@value_options)	{	if($element =~ /[<>=]/)	{	$equality="y";	}	}
				
				if($equality)
					{
					my $numopts=scalar(@value_options);
					my $numok=0;
					my $op;
					my $number;
					foreach my $el (@value_options)
						{
						if($el =~ /^([<>=]{1,2})(\s*)([+-]?[1234567890]+[$dec]?[1234567890]*e?[+-]?[1234567890]*)$/)
							{
							$op = $1;
							$number = $3;
							
							# Check if the argument value satisfies the equality/inequality
							if($op eq "<")
								{
								if($argument < $number)	{	$numok++;	}
								else	{	$message = "The argument has value $argument, but needs to be less than $number\n";	}
								}
							if($op eq ">")
								{
								if($argument > $number)	{	$numok++;	}
								else	{	$message = "The argument has value $argument, but needs to be greater than $number\n";	}
								}
							if($op eq "<=")
								{
								if($argument <= $number)	{	$numok++;	}
								else	{	$message = "The argument has value $argument, but needs to be less than or equal to $number\n";	}
								}
							if($op eq ">=")
								{
								if($argument >= $number)	{	$numok++;	}
								else	{	$message = "The argument has value $argument, but needs to be greater than or equal to $number\n";	}
								}
							
							}
						else
							{	
							die "Error in subroutine $sub_and_caller: One or more allowed values are equalities/inequalities, but the subroutine".
							" cant distinguish the operand from the number. Check formatting of the input.\n";
							}
						} # foreach element in value options (check if value satisfies equality) ends
					
					} # If there are equalities among the allowed values ends
				else
					{
					if($argument == 0)	{	$argument = "zero";	}
					unless(listTools::value_in_array($argument, $value_options_ref, "num"))	{	$message = "The argument has value $argument, but the allowed values are (@value_options)";	}
					}
				} # elseif argument is numeric ends
			} # if value options have been specified ends
		} # else (if) type is ok ends
	
	return($message);
	} # end check_arg


return(1);

# end functions