# DESCRIPTION: This script counts the number of toal and unique reads in a set of fasta files (in redundant format, NB!). It takes a textfile listing
# the files to be counted as input (one filename on each line) and writes the results to a csv outfile whose name you decide (and provide at the commandline).
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

# Define usage (error message to display if the user enters the wrong arguments at the command line)
my $usage = "Syntax error for script ${0}. Correct usage: 'perl $0 \$filenames.txt \$outfilename.csv'\n\nwhere".
"\t\$filenames.txt is a textfile listing (one on each line) the files that should be processed\n".
"\t\$outfilename.csv is the desired name of the statistics outfile\n\n";
my @pars = @ARGV or die $usage;
foreach my $el (@pars)	{	$el = text::trim($el);	}
my $namefile = shift(@pars) or die $usage;
my $outname = shift(@pars) or die $usage;
open(my $stat, ">>", $outname) or die "Script $0 couldn't create statistics outfile\n";
print($stat "File,Total,Unique\n");

# Read namefile
my @files=misc::read_list($namefile);

for(my $c=0; $c<=$#files; $c++)
	{
	my $file = $files[$c];
	my ($total, $unique)=fastaTools::count_nr_fasta_reads($file);
	print($stat "$file,$total,$unique\n");
	}

close($stat);
#close($log);
#close($wlog);
exit;

# end processing

########################################## Define local functions ##########################################

# end functions
