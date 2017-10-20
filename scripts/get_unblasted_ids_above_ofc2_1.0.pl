# DESCRIPTION: No description yet
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
my $thisfile = (__FILE__);
my $scripts = "";
if($thisfile =~ m/^(.+)\//)	{	$scripts = $1;	}
my $modules = $scripts;
$modules =~ s/scripts/modules/;

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

# Define usage (error message to display if the user enters the wrong arguments at the command line)
my $usage = "Syntax error for script ${scriptname}. Correct usage: 'perl $scriptname \$namefile \$cutoff'\n\nwhere".
"\t\$argument1 can be either 'X' or 'Y'\n".
"\t\$argument2 is a file in XX format\n".
"\t\$argument3 is the threshold value for the derived value YY to be considered credible\n\n";

my $namefile = shift or die $usage;
my $cutoff = shift or die $usage;
my $negcut = -$cutoff;
my @infiles = misc::read_filelist($namefile);
if(!open(OUT, ">>get_fasta_for_these.txt"))	{	die "Couldn't create outfile get_fasta_for_these.txt\n";	}

for(my $c=0; $c<=$#infiles; $c++)
	{
	my $genelist = $infiles[$c];
	
	# Read the relevant columns from the genelist into a matrix
	my @matrix = misc::csv_to_matrix($genelist, "all", "Feature_ID","DE_24_vs_96","OFC_24_vs_96","DE_96_vs_Lay","OFC_96_vs_Lay","type","blasted");

	my $this = scalar(@matrix);

	# Loop over matrix
	for(my $d=0; $d<=$#matrix; $d++)
		{
		my $id = $matrix[$d][0];		
		my $de1 = $matrix[$d][1];
		my $ofc1 = $matrix[$d][2];
		my $de2 = $matrix[$d][3];
		my $ofc2 = $matrix[$d][4];
		my $type = $matrix[$d][5];
		my $blasted = $matrix[$d][6];


		if($d<11)	{	print("$id\t$de1\t$ofc1\t$de2\t$ofc2\t$type\t$blasted\n");	}

		if(($de1 =~ /yes/) or ($de2 =~ /yes/))
			{
			if(($ofc1 >= $cutoff) or ($ofc1 <= $negcut) or ($ofc2 >= $cutoff) or ($ofc2 <= $negcut))
				{
				if($type =~ /novel_transcript/)
					{
					if($blasted =~ /no/)
						{
						print(OUT "$id\n");
						}
					}
				}
			}
		}
	}

close(OUT);
exit;

# end processing

########################################## Define local functions ##########################################

# end functions
