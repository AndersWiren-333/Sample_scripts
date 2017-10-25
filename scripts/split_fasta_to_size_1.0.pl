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
my $usage = "Syntax error. Usage: perl split_fasta_to_size.pl infile.fa desired_size b/kb/mb/gb\n";
my $infile = shift or die $usage;
my $desired_size = shift or die $usage;
my $unit = shift or die $usage;
if(!open(IN, $infile))	{	die "Couldn't open infile $infile\n";	}

my $desired_size_b=0;

# Check size of infile and calculate number of pieces to split it into
if($unit eq "b")	{	$desired_size_b = $desired_size;	}
elsif($unit eq "kb")        {       $desired_size_b = 1000*$desired_size;        }
elsif($unit eq "mb")        {       $desired_size_b = 1000000*$desired_size;        }
elsif($unit eq "gb")        {       $desired_size_b = 1000000000*$desired_size;        }
else	{	die "$usage;";	}

my $size_b=0;
my $part_no=1;
if(!open(OUT, ">>part_${part_no}.fa"))	{	die "Couldn't create outfile part_${part_no}.fa\n";	}

while(<IN>)
	{
	my $line = text::trim($_);
	my $temp_size = length($line);
	$size_b = $size_b+$temp_size;
	if($size_b < $desired_size_b)	{	print(OUT "$line\n");	}
	else
		{
		if($line !~ />/)	{	print(OUT "$line\n");	}
		else
			{
			close(OUT);
			$part_no++;
			$size_b=0;
			if(!open(OUT, ">>part_${part_no}.fa"))  {       die "Couldn't create outfile part_${part_no}.fa\n";     }
			print(OUT "$line\n");
			}
		}
	}

print("$part_no");

close(OUT);
close(IN);
exit;

# end processing

########################################## Define local functions ##########################################

# end functions
