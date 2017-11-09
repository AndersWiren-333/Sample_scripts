# DESCRIPTION: This script makes sure that a textfile containing one blast hit on each row is not redundant. 
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

# Declare variables and filehandles
my $usage = "Syntax error. Usage: perl de-redundize_blast_hits.pl infile.txt";
my $infile = shift or die $usage;
if(!open(IN, $infile))	{	die "Couldn't open infile $infile\n";	}
if(!open(OUT, ">>nr_$infile"))	{	die "Couldn't cretae outfile\n";	}

# Read infile into vector
my @vector=();
while(<IN>)
	{
	my $line = text::trim($_);
	push(@vector, $line);
	}

# Sort the vector
@vector = sort(@vector);

# De-redundize vector
my $ref_line=$vector[0];
print(OUT "$ref_line\n");

for(my $c=1; $c<=$#vector; $c++)
	{
	my $line = $vector[$c];
	if($line ne $ref_line)
		{
		print(OUT "$line\n");
		$ref_line = $line;
		}
	}


close(IN);
close(OUT);
exit;

# end processing

########################################## Define local functions ##########################################

# end functions
