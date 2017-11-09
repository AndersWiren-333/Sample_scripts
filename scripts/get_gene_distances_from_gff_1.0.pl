# DESCRIPTION: This script calculates distances between annotated genes in a genome, using positional information from a genome annotation file (in gff format). It estimates statistics such
# as mean, median and standard deviation for gene distance, both starnd specifically and overall. 
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

# Declare local functions (if any)
sub compute_distances;

# Define usage (error message to display if the user enters the wrong arguments at the command line)
my $usage = "Usage error for script ${0}. Correct usage: 'perl $0 \$gff_file \$exclude_negative_y_n \$offset'\n\nwhere".
"\t\$gff_file is an annotation file in gff format\n".
"\t\$exclude_negative_y_n is an indicator of whether you want to exclude negative distances. Options are 'y' and 'n'\n".
"\t\$offset is an indicator of what constitutes separation between adjacent genes (or exons) (1 nucleotide separation? 0 nucleotides?) If 0, set to 'n'\n\n";
my @pars = @ARGV or die $usage;
foreach my $el (@pars)  {       $el = text::trim($el);  }
my $gff_file = shift(@pars) or die $usage;
my $ex_neg = shift(@pars) or die $usage;
my $offset = shift(@pars) or die $usage;

if($offset eq "n")	{	$offset = 0;	}

# Open the gff file
open(my $out, ">", "gene_distances.txt") or die "Couldn't create outfile\n";

# Read the gff file into a matrix
my @matrix=fileTools::read_table($gff_file, "tsv");

# Set up new matrices to hold formatted output
my @genes_pos=();
my @genes_neg=();

my @templine=();
my $gene_counter=0;

# Loop over lines in the gff matrix
LINES: for(my $c=0; $c<=$#matrix; $c++)
	{
	# If the current line represents a gene, start a new line in the @genes matrix for that gene
	if($matrix[$c][2] eq "gene")
		{
		@templine=(); # Make sure @templine is empty
		my $gene_id="";
		if($matrix[$c][8] =~ m/(ID=gene\d+)/)	{	$gene_id = $1;	}
		else	{	print($log "Gene has no gene id\n"); next LINES;	}

		# If the start position is before the stop position
		if($matrix[$c][3] < $matrix[$c][4])
			{
			# Contig, Start, Stop, Gene_ID
			push(@templine, $matrix[$c][0], $matrix[$c][3], $matrix[$c][4], $gene_id);
			}

		# Else, if stop is before start
                elsif($matrix[$c][4] < $matrix[$c][3])
                        {
			push(@templine, $matrix[$c][0], $matrix[$c][4], $matrix[$c][3], $gene_id);
                        }
		if($matrix[$c][6] eq "+")	{	push(@genes_pos, [(@templine)]);	}
		if($matrix[$c][6] eq "-")       {       push(@genes_neg, [(@templine)]);        }
		$gene_counter++;
		}
	}
	
# Sort the @genes matrices on Contig, start, stop and Gene_ID
@genes_pos = sort { $a->[0] cmp $b->[0] || $a->[1] <=> $b->[1] || $a->[2] <=> $b->[2] || $a->[3] cmp $b->[3] } @genes_pos;
@genes_neg = sort { $a->[0] cmp $b->[0] || $a->[1] <=> $b->[1] || $a->[2] <=> $b->[2] || $a->[3] cmp $b->[3] } @genes_neg;

# Find the distances between genes
my @dists_pos=compute_distances(\@genes_pos);
my $num_neg_pos=shift(@dists_pos);
my @dists_neg=compute_distances(\@genes_neg);
my $num_neg_neg=shift(@dists_neg);

my @dists_all = (@dists_pos, @dists_neg);
my $num_neg_all = $num_neg_pos+$num_neg_neg;

# Get statistics
my $n_gene_all=scalar(@dists_all);
my $mean_gene_all=stats::mean(\@dists_all);
my $median_gene_all=stats::median(\@dists_all);
my $percentile_5_gene_all=stats::percentile(\@dists_all, 5);
my $percentile_95_gene_all=stats::percentile(\@dists_all, 95);
my $range_gene_all=stats::range_string(\@dists_all);
my $stdev_gene_all=stats::stdev(\@dists_all);
my $two_stdev_gene_all = 2*$stdev_gene_all;

my $n_gene_pos=scalar(@dists_pos);
my $mean_gene_pos=stats::mean(\@dists_pos);
my $median_gene_pos=stats::median(\@dists_pos);
my $percentile_5_gene_pos=stats::percentile(\@dists_pos, 5);
my $percentile_95_gene_pos=stats::percentile(\@dists_pos, 95);
my $range_gene_pos=stats::range_string(\@dists_pos);
my $stdev_gene_pos=stats::stdev(\@dists_pos);
my $two_stdev_gene_pos = 2*$stdev_gene_pos;

my $n_gene_neg=scalar(@dists_neg);
my $mean_gene_neg=stats::mean(\@dists_neg);
my $median_gene_neg=stats::median(\@dists_neg);
my $percentile_5_gene_neg=stats::percentile(\@dists_neg, 5);
my $percentile_95_gene_neg=stats::percentile(\@dists_neg, 95);
my $range_gene_neg=stats::range_string(\@dists_neg);
my $stdev_gene_neg=stats::stdev(\@dists_neg);
my $two_stdev_gene_neg = 2*$stdev_gene_neg;


# Print statistics to an outfile
print($out "Gene-Gene distances (all genes) (n = $n_gene_all):\n");
if($ex_neg eq "y")	{	print($out "$num_neg_all negative distances (excluded)\n");	}
else	{       print($out "$num_neg_all negative distances (included)\n");     }
print($out "Mean: $mean_gene_all\n");
print($out "Median: $median_gene_all\n");
print($out "5th percentile: $percentile_5_gene_all\n");
print($out "95th percentile: $percentile_95_gene_all\n");
print($out "Range: $range_gene_all\n");
print($out "Standard deviation: $stdev_gene_all\n");
print($out "2*Standard deviation: $two_stdev_gene_all\n");
print($out "\n\n");

print($out "Gene-Gene distances (positive strand) (n = $n_gene_pos):\n");
if($ex_neg eq "y")     {       print($out "$num_neg_pos negative distances (excluded)\n");     }
else    {       print($out "$num_neg_pos negative distances (included)\n");     }
print($out "Mean: $mean_gene_pos\n");
print($out "Median: $median_gene_pos\n");
print($out "5th percentile: $percentile_5_gene_pos\n");
print($out "95th percentile: $percentile_95_gene_pos\n");
print($out "Range: $range_gene_pos\n");
print($out "Standard deviation: $stdev_gene_pos\n");
print($out "2*Standard deviation: $two_stdev_gene_pos\n");
print($out "\n\n");

print($out "Gene-Gene distances (negative strand) (n = $n_gene_neg):\n");
if($ex_neg eq "y")     {       print($out "$num_neg_neg negative distances (excluded)\n");     }
else    {       print($out "$num_neg_neg negative distances (included)\n");     }
print($out "Mean: $mean_gene_neg\n");
print($out "Median: $median_gene_neg\n");
print($out "5th percentile: $percentile_5_gene_neg\n");
print($out "95th percentile: $percentile_95_gene_neg\n");
print($out "Range: $range_gene_neg\n");
print($out "Standard deviation: $stdev_gene_neg\n");
print($out "2*Standard deviation: $two_stdev_gene_neg\n");
print($out "\n\n");




close($out);
close($log);
#close($wlog);
exit;

# end processing

########################################## Define local functions ##########################################

sub compute_distances
	{
	my $matrixref = $_[0];
	my @matrix = @{$matrixref};
	my $num_neg=0;

	my @gene_distances=();

	my $comp_contig=$matrix[0][0];
	my $comp_gene_start="banana";	# Setting this to an arbitrary string will give us an error message later if the value fails to be updated as it should below (otherwise we will never know if we have used the worong values)
	my $comp_gene_stop="orange";

	# Loop over genes
	for(my $cc=0; $cc<=$#matrix; $cc++)
		{
		my $contig = $matrix[$cc][0];

		my $gene_start=$matrix[$cc][1];
		my $gene_stop=$matrix[$cc][2];
		my $last = (scalar(@{$matrix[$cc]}))-1;

		# If this is the first gene or id, then we are on a new contig. Don't compute intergene distance
		if(($cc==0) or ($contig ne $comp_contig))
			{
			$comp_contig=$contig;
			$comp_gene_start=$matrix[$cc][1];
			$comp_gene_stop=$matrix[$cc][2];
		
			} # If first gene ends

		# If this is not the first gene in on the list or if we are still on the same contig
		else
			{
			# Compute gene-gene distance
			my $gene_start = $matrix[$cc][1];
			my $gene_stop = $matrix[$cc][2];

			my $gene_distance = $gene_start-$comp_gene_stop-$offset;
			if($gene_distance < 0)	{	$num_neg++;	}
                       	if($ex_neg eq "y")
                        	{
                               	if($gene_distance>0)    {       push(@gene_distances, $gene_distance);  }
                              	}
			else	{	push(@gene_distances, $gene_distance);	}
			($comp_gene_start, $comp_gene_stop) = ($gene_start, $gene_stop);

			} # If not first gene ends
		} # Loop over genes ends

	unshift(@gene_distances, $num_neg);
	return(@gene_distances);
	}


# end functions
