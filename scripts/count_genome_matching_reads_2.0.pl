# DESCRIPTION: This perl script counts the number of redundant and non-redundant reads in a patman output file (genome alignment results, "filename.pat")
# and prints the results on the screen (the results can also be captured and processed if this script is called from another script).
# end description

########################################## Universal perl script header ##########################################

#!/usr/bin/perl

# perl_script_update

# Load libraries that this script depends on
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
my $maintain = $FindBin::Bin;
my $scripts = $maintain;
$scripts =~ s/\/maintainance/\/scripts/;
my $modules = $scripts;
$modules =~ s/\/scripts/\/modules/;

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
require "$modules/combinatorics.pm";
require "$modules/db.pm";

# Create a timestamp string (can be attached to the name of logfiles, for example
my $timestamp = envir::timestamp();
my $rscript = "Rscript";

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

# Declare variables and filehandles
my $usage = "Syntax error. Usage: perl count_genome_matching_reads.pl infile.pat";
my $infile = shift or die $usage;

if(!open(IN, "$infile"))	{	die "Couldn't open $infile";	}
my @arr=();

while(<IN>)
	{
	my $line = text::trim($_);
	if($line =~ /([ACTG]+-\d+)/)
		{
		my $seq = $1;
		$seq = text::trim($seq);
		push(@arr, $seq);
		}
	}


my @arrS = sort(@arr);
my $comp = $arrS[0];

my $countU=0;
my $countR=0;

if(defined($comp))
	{
	$countU = 1;
	#my $countR = 0;
	if($comp =~ /-(\d+)/)	{	$countR = $1;	}

	# Loop over elements in @arrS
	for(my $i=0 ; $i<=$#arrS; $i++)
		{
		if( $arrS[$i] eq $comp )	{	next;	}	# Skip identical reads (the same read may align to several places in the genome, but should only be counted once).
		else	
			{
			$countU++;
			if($arrS[$i] =~ /-(\d+)/)
				{
				my $num=$1;
				$countR=$countR+$num;
				}
			$comp=$arrS[$i]; 
			}
		}
	
	}
print("$countR,$countU\n");

exit;

# end processing

########################################## Define local functions ##########################################

# end functions
