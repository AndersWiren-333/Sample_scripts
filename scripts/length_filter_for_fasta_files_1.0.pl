# DESCRIPTION: This is a template for writing new perl scripts
# end description

########################################## Universal perl script header ##########################################

#!/usr/bin/perl

# Load libraries that this script depends on
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

# Load home-made modules (aka "libraries", aka "packages") that this script depends on
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

# Create a logfile
#open(my $log, ">>", "log_${0}_${timestamp}.txt") or die "Script $0 couldn't create logfile\n";

# Create warninglog
#open(my $wlog, ">>", "warninglog_${0}_${timestamp}.txt") or die "Script $0 couldn't create warninglog\n";
#local $SIG{__WARN__} = sub
#       {
#       my $message = shift;
#       print($wlog $message);
#       };

# end header

########################################## Processing ##########################################

# Declare local functions (if any)

# Define usage (error message to display if the user enters the wrong arguments at the command line)
my $usage = "Usage error for script ${0}. Correct usage: 'perl $0 \$infile.fa \$length_cutoff, \$outname'\n\nwhere".
"\t\$infile.fa is the file that should be used as input\n".
"\t\$length_cutoff is the minimum number of basepairs a sequence needs to be to be retained in the outfile\n".
"\t\$outname is the desired name of the outfile that will be produced\n\n";
my @pars = @ARGV or die $usage;
foreach my $el (@pars)  {       $el = text::trim($el);  }
my $infile = shift(@pars) or die $usage;
my $length_cutoff = shift(@pars) or die $usage;
my $outname = shift(@pars) or die $usage;

# Read the file
my @matrix=fastaTools::fasta_to_matrix($infile);

# Loop over reads in the file and exclude all that are shorter than $length_cutoff
open(my $out, ">", $outname) or die "Script $0 couldn't create outfile $outname\n";

for(my $c=0; $c<=$#matrix; $c++)
	{
	if(length($matrix[$c][1]) >= $length_cutoff)
		{
		print($out ">$matrix[$c][0]\n$matrix[$c][1]\n");
		}
	}

close($out);
#close($log);
#close($wlog);
exit;

# end processing

########################################## Define local functions ##########################################

# end functions
