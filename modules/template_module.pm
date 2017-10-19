package template_module;

# sub sub_1($parameter)
# sub sub_2(\@matrixref)
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
my $scripts = cwd;
chdir("../modules");
my $modules = cwd;
chdir("../scripts");

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

# Description of what the function does and how to use it
sub sub_1
	{
	# Set error messages and accept input parameters
        my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
        my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$argument1, \$argument2, \@arguments3)'\n\nwhere".
        "\t\$argument1 can be either 'X' or 'Y'\n".
        "\t\$argument2 is a file in XX format\n".
        "\t\@arguments3 is a list of arguments\n\n";
        my @pars = @_ or die $usage;
        foreach my $el (@pars)  {       $el = text::trim($el);  }
        my $argument1 = shift @pars or die $usage;	# NB! This will not work if the argument is the number 0 (because it will then be interpreted as false). In that case you need to use 'shift(@pars) or 0', but it may lead to problems....
        my @other_arguments = @pars or die $usage;

	open(my $in, "<", $infile) or die "Subroutine $subname (called by script '${calling_script}', line ${calling_line}) couldn't open infile $infile\n";
	open(my $out, ">", $outname) or die "Subroutine $subname (called by script '${calling_script}', line ${calling_line}) couldn't create outfile $outname\n";

	close($in);
	close($out);
	return($something);
	}
	
return(1);

# end functions
