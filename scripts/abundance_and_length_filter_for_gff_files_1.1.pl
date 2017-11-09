# DESCRIPTION: This script takes as input a tab separated list of information about novel (unannotated) transcripts (transcripts assembled from an RNAseq data set that haven't been annotated
# before in the species of interest) (= an raw_una.gff file). It filters out those transcripts that are shorter than a certain length and have less than a certain abundance in the
# dataset as a whole (these values are given as arguments at the command line). It writes all transcripts that passed the filter to a tab separated outfile (= una.gff) with the
# same format as the input file. 
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
#use diagnostics;
use FindBin;
use DBI;
use IO::Compress::Gzip qw(gzip);
use IO::Uncompress::Gunzip qw(gunzip);

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

# Define usage (error message to display if the user enters the wrong arguments at the command line)
my $usage = "Syntax error for script ${0}. Correct usage: 'perl $0 \$infile.gff \$length_cutoff \$abundance_cutoff \$outfile_name'\n\nwhere".
"\t\$infile.gff is a file listing chromosome, start and stop position of features\n".
"\t\$length_cutoff is the minimum length in basepairs a feuture needs to have for inclusion in the filtered outfile\n".
"\t\$abundance_cutoff is the minimum abundance a feature needs to have to be included in the filtered outfile\n".
"\t\$outfile_name is the desired name of the outfile\n\n";

# Declare variables and filehandles
my $infile = shift or die $usage;
my $length_criterion = shift or die $usage;
my $abundance_criterion = shift or die $usage;
my $outname = shift or die $usage;

open(my $in, "<", $infile) or die "Script $0 couldn't open infile $infile\n";
open(my $out, ">>", $outname) or die "Script $0 couldn't create outfile $outname\n";

# Format of a line in the infile:
#	0	seqid				NC_015762.1
#	1	start				311
#	2	stop				2062
#	3	total abundance			551
#	4	positive strand abundance	548
#	5	negative strand abundance	3

# Loop over lines in the infile
while(my $line = <$in>)
	{
	my $line = text::trim($line);
	my @arr=split("\t", $line);
	foreach my $el (@arr)	{	$el = text::trim($el);	}
	my $length = $arr[2]-$arr[1]+1;
	if($length >= $length_criterion && $arr[3] >= $abundance_criterion)	{	print($out "$line\n");	}
	}

close($in);
close($out);

#close($log);
#close($wlog);
exit;

# end processing

########################################## Define local functions ##########################################

# end functions
