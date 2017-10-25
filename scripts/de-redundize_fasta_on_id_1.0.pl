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

# Define usage (error message to display if the user enters the wrong arguments at the command line)
my $usage = "Syntax error for script ${0}. Correct usage: 'perl $0 \$infile1 \$infile2 etc.'\n\nwhere".
"\t\$infile1 is the first file that should be used as input\n".
"\t\$infile2 is the second file that should be used as input\n\n";
my @files = @ARGV or die $usage;
my @all_matrix=();

for(my $c=0; $c<=$#files; $c++)
	{
	my $file = $files[$c];
	my @matrix=fastaTools::fasta_to_matrix($file);
	push(@all_matrix, @matrix);
	@matrix=();
	}

my @unique_matrix=misc::unique_matrix(\@all_matrix, 1, "alph", "n");
@all_matrix=();

if(!open(OUT, ">>nr_refseqs.fa"))	{	die "$scriptname couldn't create outfile nr_refseqs.fa\n";	}
for(my $d=0; $d<=$#unique_matrix; $d++)
	{
	my $id=$unique_matrix[$d][0];
	my $seq=$unique_matrix[$d][1];	
	print(OUT ">$id\n$seq\n");
	}
close(OUT);
exit;

# end processing

########################################## Define local functions ##########################################

# end functions
