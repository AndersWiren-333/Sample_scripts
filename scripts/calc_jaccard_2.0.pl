# DESCRIPTION: Calculates jaccard coefficients
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

# Open filehandles, declare variables
my $usage = "perl calc_jaccard.pl gene_expression_matrix.csv limit sample_name_1 sample_name_2 etc...";
my $infile = shift or die $usage;
my $limit = shift or die $usage;
my @names=("Sample", @ARGV);

if(!open(IN, $infile))	{	die "Couldn't open infile";	}
if(!open(OUT, ">>jaccard_top_${limit}_genes.csv"))    {       die "Couldn't create outfile";  }

my $namn=join(",", @names);
my $header=("$namn\n");
print(OUT $header);

# Read raw data from infile into a 2-dimensional array

my @matrix=();
my $cols=0;

while(<IN>)
	{
	my $line = $_;
	my @arr=split(",", $line);
	for(my $j=0; $j<=$#arr; $j++)	{	$arr[$j] = text::trim($arr[$j]);	}
	push(@matrix, \@arr);
	$cols=scalar(@arr);
	}

# Find the top X most abundant genes for each sample

for(my $i=1; $i<$cols; $i++)	# For each data column, do this
	{
	@matrix = sort {$a->[$i] <=> $b->[$i]} @matrix;
	my $rows=scalar(@matrix);
	for(my $rad=0; $rad<$rows; $rad++)
		{	
		if($rad<$limit)	{	$matrix[$rad][$i]=1;	}	# If we are below the limit (e.g. top 500 genes)... set the value for the current gene in the crrent sample to 1(= it's present)
		else	{	$matrix[$rad][$i]=0;	}		# Set all other genes to 0 (= not present on top 500 list)
		}
	}

# Compare each sample (column) in the array to each other and calculate the Jaccard index between them

for(my $i=1; $i<$cols; $i++)	# For each column ($i, sample), do this...
	{
	my @line=();
	
	# Compare to each other sample ($j)...
	for(my $j=1; $j<$cols; $j++)
		{
		#print("\nComparing $names[$i] to $names[$j] ...");
		my $intersection=0;
		my $union=0;
		
		# For each line (gene), determine if that line contributes to the intersection and union
		foreach my $line (@matrix)
			{
			if( ($line->[$i])!=0 and ($line->[$j])!=0 )     {       $intersection++;     }
			if( ($line->[$i])!=0 or ($line->[$j])!=0)      {       $union++;    }
			}
		my $jaccard=($intersection/$union);
		push(@line, $jaccard);
		}
	my $line2=join(",", @line);
	my $outline=("$names[$i],".$line2."\n");	
	print(OUT $outline);
	}

print("\n\nFinished!\n");
close(IN);
close(OUT);

exit;

# end processing

########################################## Define local functions ##########################################

# end functions
