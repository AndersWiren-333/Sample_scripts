# DESCRIPTION: This script takes a *.txt outfile from FastQC and extracts information about the basecall quality from it. The script then
# computes the average basecall quality score over all bases and reads in the original *.fastq file.
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
use FindBin;

# Set paths to scripts and modules. Setting explicit paths to the scripts and modules in this specific repository (rather than adding paths to @INC, PERLLIB and PATH on
# your local system) avoids the risk of scripts calling the wrong scripts/modules if you have other repositories on your system that happen to have some script- and module names
# in common with this repository.
my $scripts = $FindBin::Bin;
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

# Declare variables and filehandles
my $usage = "perl extract_quality.pl infile.txt";
my $infile = shift or die $usage;
if(!open(IN, $infile))	{	die "Couldn't open infile $infile\n";	}

my @cont = <IN>;			# Read all lines of the infile into an array
my $content = join("", @cont);		# And make sure you have theentire contents as one string

my $total_reads=0;
my $avg_qual=0;

# Get total number of reads
if($content =~ /(Total Sequences\s+)([\d]+)(Sequences flagged as poor quality)/)	{	my $total_reads = $2;	}

if($content =~ /(>>Per base sequence quality\s+pass\n)([\s\w\W\d\D]+)(>>END_MODULE\n>>Per tile sequence quality)/)
	{
	my $ex = $2;
	my @big_arr = split("\n", $ex);
	foreach my $el (@big_arr)	{	$el = text::trim($el);	}
	shift(@big_arr);
	my $totscore=0;
	my $count=0;

	for(my $c=0; $c<=$#big_arr; $c++)
		{
		my $line = $big_arr[$c];
		$line =~ s/\s/H/g;
		$line =~ s/H+/,/g;
		my @arr = split(",", $line);
		foreach my $el (@arr)	{	$el = text::trim($el);	}
		

############################################ 	HAR AR JAG NU! ##################################################

		}

	foreach my $line (@crudearr)
		{
		$count++;
		if($line =~ /(\d+)(\s+)(\d+.\d+)(\s+)(.+)/)
			{
			my $qscore = $3;
			$totscore = ($totscore + $qscore);
			}		
		}
	$avg_qual=($totscore/$count);
	}

my $prob_incorr=10^(-0.1*$avg_qual)


print("$avg_qual,$prob_incorr,$total_reads");

close(IN);
exit;

# end processing

########################################## Define local functions ##########################################

# end functions
