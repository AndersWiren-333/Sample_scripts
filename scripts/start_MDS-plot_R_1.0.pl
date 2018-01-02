# DESCRIPTION: This script takes a gene expression matrix (in csv format) and adds a header containing sample names to it. It calculates the average
# expression for each gene across samples (across columns) and adds that value as an additional column. It then sorts the genes from highest
# to lowest average expression, after which it uses the MDS_plot_2.0.R to produce multidimensional scaling (MDS) plots for the data.
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

# Declare variables and filehandles
my $usage="Usage error for script ${0}. Correct usage: 'perl $0 \$expression_matrix.csv \$sample_name_01 \$sample_name_02 ... \$colour_sample_1 \$colour_sample_2 etc...'\n\nwhere".
"\t\$expression_matrix.csv is a gene expression matrix file - without headers - in csv format\n".
"\t\$sample_name_01, \$sample_name_02 etc. are the names of the samples/libraries in the dataset\n".
"\t\$colour_sample_01, \$colour_sample_02 etc. are the colours that samples should be represented by in the output plots\n\n";
my $infile = shift or die $usage;
if(!open(IN, $infile))	{	die "Couldn't open $infile";	}

my @names_and_cols=@ARGV;
my $num=scalar(@names_and_cols);
my $fetch=$num/2;
my @names=(@names_and_cols[0..($fetch-1)]);
my @colours=(@names_and_cols[$fetch..$#names_and_cols]);

my @header=("Gene", @names);

my $out_header=join(",", @header);
my $last_print=$#header;
my $mean_index=($last_print)+1;
my @matrix=();

# Read in the expression matrix line by line, calculate average expression as you go along

while(<IN>)
	{
	my $line = $_;
	my @arr=split(",", $line);
	for(my $a=0; $a<=$#arr; $a++)	{	$arr[$a] = text::trim($arr[$a]);	}

	my $total=0;
	my $n=(scalar(@arr))-1;			# The sequence id constitutes one value and shouldn't be counted, hence the "-1"
	for(my $i=1; $i<scalar(@arr); $i++)
		{
		$total=($total)+($arr[$i]);		
		}
	my $mean=($total)/($n);
	push(@arr, $mean);
	push(@matrix, \@arr);
	}

close(IN);

# Sort the matrix (descending) on average expression
my @s_matrix = sort { $b->[$mean_index] <=> $a->[$mean_index] } @matrix;

my $num_genes=scalar(@s_matrix);
my $k=(int($num_genes/1000))+1;
print("\n\n\tNumber of parts: $k\n\n");
my $j=0;


# Print out csv files containing the top most expressed 1-1000, 1001-2000 etc. genes
for(my $part=1; $part<=$k; $part++)
        {
        if(!open(OUT, ">>$part.csv"))   {       die "Couldn't create outfile";  }
        print(OUT "$out_header\n");
        my $lim=1000*$part;
        while($j<$lim)
                {
                if($j==scalar(@s_matrix))    {       last;   }	# If this is the line after the last line of @s_matrix, don't do anything more.
		my $arref=$s_matrix[$j];
		my @arr=@{$arref};
		my $outline=join(",", @arr[0..$last_print]);
                print(OUT "$outline\n");
                $j++;
                }
        close(OUT);
        }

# Copy the infile and name the copy for easy use in script scale.R
if(!open(ALL, ">>all.csv"))	{	die "Couldn't create all.csv\n";	}
print(ALL "$out_header\n");

for(my $a=0; $a<=$#s_matrix; $a++)
	{
	my @arr=@{$s_matrix[$a]};
	my $outline=join(",", @arr[0..$last_print]);
	print(ALL "$outline\n");
	}
close(ALL);

# Run scale.R
system("Rscript $scripts/MDS_plot_2.0.R $k $infile @colours");

# Remove intermediate csv files
for(my $b=1; $b<=$k; $b++)
	{
	my $file=("$b.csv");
	unlink($file);
	}
unlink("all.csv");
exit;

# end processing

########################################## Define local functions ##########################################

# end functions
