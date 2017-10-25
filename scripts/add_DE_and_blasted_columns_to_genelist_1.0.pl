# DESCRIPTION: This script adds extra columns to a list of differentially expressed genes; one column stating whether a gene/feature is differentially expressed or not ("yes"/"no",
# one column for each comparison being made) and one column stating whether or not the gene/feature sequence has been blasted to search for homologous sequences outside
# the organism of interest ("yes"/"no"). The script takes as input the genelist in question and a list of feature IDs of the features that have been blasted (text file, one ID on each row, ID format
# "NC_1234.5_678_91011").
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
my $usage = "perl add_DE_and_blasted_columns_to_genelist.pl genelist.csv blasted_genes.txt\n";
my $genelist = shift or die $usage;
my $blastlist = shift or die $usage;
if(!open(GENES, $genelist))	{	die "Couldn't open genelist $genelist\n";	}
if(!open(BLAST, $blastlist))	{	die "Couldn't open blastlist $blastlist\n";	}

# Read blastlist into vector
my @blasted=();

while(<BLAST>)
	{
	my $line = trim($_);
	push(@blasted, $line);
	}
close(BLAST);
@blasted = sort { $a cmp $b } @blasted;


# Read genelist into matrix
if(!open(OUT, ">>new_$genelist"))	{	die "Couldn't create outfile new_$genelist\n";	}

my @headers=();
my @new_headers=();
my $genecounter=0;
my @ofc_indices=();
my $blast_index_1=0;
my $lastOFCplusone=0;

while(<GENES>)
	{
	# read current line into memory, split into columns and make sure there is no hidden whitespace
	my $line = trim($_);
	my @arr = split(",", $line);
	foreach my $el (@arr)	{	$el = trim($el);	}

	# Extract headers
	if($genecounter==0)
		{
		@headers = @arr;
		(@ofc_indices) = grep { $headers[$_] =~ /OFC/ } 0..$#headers;
		($blast_index_1) = grep { $headers[$_] eq "BLAST_hit_01" } 0..$#headers;
		push(@new_headers, @headers[0..(($ofc_indices[0])-1)]);				# Add columns Feature_ID, raw data and length to new headers

		# For each comparison (i.e. for each ofc column)
		for(my $c=0; $c<=$#ofc_indices; $c++)
			{
			my $index = $ofc_indices[$c];	# This is the position of the ofc column in @headers
			my $name = $headers[$index];	# This is the name of the ofc column
			$name =~ s/OFC/DE/;		# Make the name of a new column
			push(@new_headers, $name, $headers[$index]);		# Add the name of the new column and the name of the ofc column to @new_headers
			}
		$lastOFCplusone = (scalar(@new_headers))-2;
	
		my $mid_end = ($blast_index_1)-1;
		my @header_mid = @headers[$lastOFCplusone..$mid_end];

		push(@new_headers, @header_mid, "blasted", @headers[$blast_index_1..$#headers]);
		my $outheaders = join(",", @new_headers);
		print(OUT "$outheaders\n");
		}
	
	# Process line
	else
		{
		my $seqid = $arr[0];
		my @new_arr = @arr[0..(($ofc_indices[0])-1)];
		for(my $d=0; $d<=$#ofc_indices; $d++)
                        {
                        my $index = $ofc_indices[$d];
                        my $ofc = $arr[$index];
			my $de = "no";
			if($ofc != 0)	{	$de = "yes";	}
                        push(@new_arr, $de, $arr[$index]);
                        }
		push(@new_arr, @arr[$lastOFCplusone..(($blast_index_1)-1)]);

		# Loop through list of blasted features to see if current feature has been blasted
		my $is_blasted="no";

		for(my $e=0; $e<=$#blasted; $e++)
			{
			my $blast_id = $blasted[$e];
			if($seqid eq $blast_id)	{	$is_blasted = "yes";	}
			}

		push(@new_arr, $is_blasted, @arr[$blast_index_1..$#arr] );
		
		my $outline = join(",", @new_arr);
		print(OUT "$outline\n");
		}
	$genecounter++;
	}

close(GENES);
close(OUT);
exit;

# end processing

########################################## Define local functions ##########################################

# end functions
