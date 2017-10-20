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

# Declare variables and filehandles
my $usage = "perl generate_novel_reads.pl transcript_length read_length num_transcripts num_reads";
my $tra_length = shift or die $usage;
my $read_length = shift or die $usage;
my $num_transcripts = shift or die $usage;
my $num_reads = shift or die $usage;
if(!open(READS, ">>new_reads.fastq"))	{	die "Couldn't create outfile\n";	}
if(!open(TRANS, ">>new_transcripts.fa"))	{	die "Couldn't create another outfile\n";	}

my $nu=1;

# Loop over number of transcripts
for(my $f=1; $f<=$num_transcripts; $f++)
	{
	my @transcript=();
	# Build a transcript (Loop over number of nucleotides)
	for(my $c=1; $c<=$tra_length; $c++)
		{
		my $number=rand(10);		
		if($number<=2.5)	{	push(@transcript, "A");	}
		elsif($number<=5 && $number>2.5)        {       push(@transcript, "C");	}
		elsif($number<=7.5 && $number>5)        {       push(@transcript, "G");	}
		elsif($number<=10 && $number>7.5)        {       push(@transcript, "T");	}
		}

	my $trans=join("", @transcript);
	print(TRANS ">AW_111111.${nu} Bombus rex\n$trans\n");
	$nu++;
	# Break transcript into x overlapping pieces of length $read_length
	my @reads=();
	my $start=0;
	for(my $d=1; $d<=$num_reads; $d++)
		{
		my $read=substr($trans, $start, $read_length);
		push(@reads, $read);
		$start=$start+($read_length-5);
		if($tra_length-$start<$read_length)	{	$start=0;	}
		}

	my $qual_line="J" x $read_length;

	# Print reads
	for(my $e=0; $e<=$#reads; $e++)
		{
		print(READS "\@HWI-D00148:546:H77GTBCXY:1:2213:12952:79733/1\n");
		print(READS "$reads[$e]\n");
		print(READS "+\n");
		print(READS "$qual_line\n");
		}
	}

close(READS);
close(TRANS);
exit;

# end processing

########################################## Define local functions ##########################################

# end functions
