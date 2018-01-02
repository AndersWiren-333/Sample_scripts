# DESCRIPTION: This script takes a sequence read file in fastq format and converts the quality scores of each sequence to numbers,
# which are then written to a csv outfile (only the quality lines of the fastq infile).
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
# use diagnostics;

# Set path to custom modules
unshift(@INC,
"/scratch2/Anders/01_Projects/01_Scripts/00_Modules", "E:\\UEA\\01_Projects\\01_Scripts\\00_Modules",
"/Users/anders2/Documents/03_UEA/01_Scripts/00_Modules", "/gpfs/home/zjh15phu/01_Scripts/00_Modules");

# Load home-made modules (aka "libraries", aka "packages") that this script depends on
require envir;
require dnaTools;
require fastaTools;
require fastqTools;
require matrixTools;
require misc;
require rtf;
require stats;
require text;

# Get environment information
my ($timestamp, $scripts, $r, $rscript, $compress, $uncompress) = envir::senseEnv();

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

# Declare variables and file handles
my $usage  = "perl quality_plot.pl infile.fastq";
my $infile = shift or die $usage;
if(!open(IN, $infile))	{	die "Couldn't open $infile\n";	}
if(!open(OUT, ">qplot_$infile.csv"))	{	die "Couldn't create outfile\n";	}

my $offset = 33;

# Define quality scale and limit
#my $scale="!\"#$%&'()*+,-./0123456789:;<=>?\@ABCDEFGHIJKLMNOPQRSTUVWXYZ";
#my @symbols=split("", $scale);
#my @numbers=(0..56);
#print("Symbols: ".scalar(@symbols)."\tNumbers: ".scalar(@numbers)."\n");

my $line_counter = 0;

# Loop over sequences
while(<IN>)
	{
	$line_counter++;
	if($line_counter%4 == 0)        
		{
		my $line = text::trim($_);

		for($i = 0; $i < length($line)-1; $i++)
			{
			$ascii_char = (ord(substr($line,$i,1)))-$offset;
			print OUT "$ascii_char,"; 
			}
		$ascii_char = (ord(substr($line,$i,1)))-$offset;
		print OUT "$ascii_char\n";
		}
	}

close(IN);
close(OUT);
exit;

# end processing

########################################## Define local functions ##########################################

# end functions
