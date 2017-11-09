# DESCRIPTION: This script takes an expression matrix (a table of read counts for a set of features/genes/reads, one on each row, in a number of samples/libraries,
# one in each column) and filters out short and lowly expressed features. The exclusion is based on expression within groups of replicate samples/libraries that are
# defined in a sample name textfile. Only those features/genes/reads will be retained where
#
#	a) at least a specific number ($w_groups) of replicate sample groups (minimum is 1)
#	b) have at least a specific number ($x_samples) or percentage ($x_perc_samples) of the samples (if $x_samples samples should be used, set $x_perc_samples to "n" and vice-versa)
#	c) where the read count is greater than or equal to a certain threshold ($y_readcount)
#
# An optional additional criterion for a feature to be retained is that the same sample group that satisfies the other criteria also should have an average read
# count of at least a specified number ($z_expr_mean). If this criterion should not be used, set that parameter to "n", otherwise set it to a number.  
#
# The sample name text file must list one sample on each line, the sample name being immediately followed by a comma and the name (NB! Not number!) of the group that
# sample belongs to, e.g:
#
#		sample_1,worker
#		sample_2,worker
#		sample_3,queen
#		sample_4,queen
#
# Example usage:
# perl abundance_and_length_filter_for_matrices.pl gene_expression_matrix_2.csv sample_names.txt filtered_expression_matrix.csv 19 1 2 "n" 20 "n"
#
# (filter "gene_expression_matrix_2.csv" using sample name information from "sample_names.txt", let the filtered outfile be called "filtered_expression_matrix.csv". Retain only features that are at
# least 19 basepairs long and for which at least 1 sample group has at least 2 samples that have a read count of at least 20. Don't use the "$x_perc_samples" or "$z_expr_mean" criteria.)
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
my $usage1 = "Usage error for script ${0}.\n\n";

my $usage2 = "Correct usage: 'perl $0 \$gene_expression_matrix.csv \$sample_names.txt \$outfile_name \$length_limit \$w_groups \$x_samples \$x_perc_samples \$y_readcount \$z_expr_mean\n\nwhere".
"\t\$gene_expression_matrix.csv is the expression matrix to be filtered\n".
"\t\$sample_names.txt is a plain text file with, on each line, the name of a sample followed by comma and the the name of the group that sample belongs to\n".
"\t\$outfile_name is the desired name for the filtered outfile\n".
"\t\$length_limit is the minimum length (in bp) a feature should have to be retained\n".
"\t\$w_groups is the number of sample groups that need to satisfy the read count criteria for the feature to be retained.\n".
"\t\$x_samples is the number of samples in each of \$w_groups sample groups that need to have a read count above a specific threshold. If percentage of samples should be used instead, set this option to 'n'.\n".
"\t\$x_perc_samples is the percentage of samples in each of \$w_groups sample groups that need to have a read count above a specific threshold. If number of samples should be used instead, set this option to 'n'.\n".
"\t\$y_readcount is the readcount that \$x_samples need to be greater than or equal to for a feature to be retained.\n".
"\t\$z_expr_mean is the minimum mean readcount that the samples in each of \$w_groups sample groups need to have for a feature to be retained\n\n";

my $usage3 = "Did you try to give 0 (zero) as a value for a numerical argument to this script (Perl can't handle that)? Or did you forget to supply a value?\n\n";

my $usage = $usage1.$usage2;
my $usage_numer = $usage1.$usage3.$usage2;

my @pars = @ARGV or die $usage;
foreach my $el (@pars)  {       $el = text::trim($el);  }
my $expression_matrix = shift(@pars) or die $usage;
my $sample_names_file = shift(@pars) or die $usage;
my $outname = shift(@pars) or die $usage;
my $length_cutoff = shift(@pars) or die 0;
my $w_groups_cutoff = shift(@pars) or 0;
my $x_samples_cutoff = shift(@pars) or 0;
my $x_perc_samples_cutoff = shift(@pars) or 0;
my $y_readcount_cutoff = shift(@pars) or 0;
my $z_expr_mean_cutoff = shift(@pars) or 0;

open(my $in, "<", $expression_matrix) or die "Script $0 couldn't open infile $expression_matrix\n";
open(my $out, ">", $outname) or die "Script $0 couldn't create outfile $outname\n";

# Create an outname
my $proposed_outname="";
if($x_perc_samples_cutoff eq "n")
	{
	if($z_expr_mean_cutoff eq "n")	{	$proposed_outname = "filt_le${length_cutoff}_gr${w_groups_cutoff}_s${x_samples_cutoff}_ab${y_readcount_cutoff}_${expression_matrix}";	}
	else	{	$proposed_outname = "filt_le${length_cutoff}_gr${w_groups_cutoff}_s${x_samples_cutoff}_ab${y_readcount_cutoff}_mab${z_expr_mean_cutoff}_${expression_matrix}";	}
	}
elsif($x_samples_cutoff eq "n")
	{
        if($z_expr_mean_cutoff eq "n")  {       $proposed_outname = "filt_le${length_cutoff}_gr${w_groups_cutoff}_spc${x_perc_samples_cutoff}_ab${y_readcount_cutoff}_${expression_matrix}";   }
        else    {       $proposed_outname = "filt_le${length_cutoff}_gr${w_groups_cutoff}_spc${x_perc_samples_cutoff}_ab${y_readcount_cutoff}_mab${z_expr_mean_cutoff}_${expression_matrix}";  }
	}

# Read the sample names file into a matrix of group indices
my @group_indices=text::read_sample_list($sample_names_file);

# Read in the expression matrix infile, line by line, while simultaneously filtering on abundance and length
while(my $line = <$in>)
	{
	# Make sure the current line and its elements don't contain hidden whitespace that may lead to erroneous interpretation of
	# numbers and textstrings
	my $line = text::trim($line);
	my @arr=split(",", $line);
	foreach my $el (@arr)	{	$el=text::trim($el);	}
	
	# Calculate the length of the feature - either from information in the feature name/identifier (if the Feature_ID format is 'NC_123.4_startpos_stoppos')
	# or from the length of the feature ID (if the Feature_ID format is 'ACTGGCTA...')	
	my $length=0;
	if($arr[0] =~ /^[ACTGU]+$/m)	{	$length = length($arr[0]);	}
	else
		{
		my @feat_name = split("_", $arr[0]);
		$length=($feat_name[3])-($feat_name[2])+1;
		}

	# Check if the feature passes the length criterion
	my $pass_length="no";
	if($length >= $length_cutoff)	{	$pass_length="yes";	}
	my @data = @arr;
	shift(@data);

	# Check how many of the sample groups pass the read count criteria
	my $pass_abund_criteria=0;

	# Loop over groups
	for(my $c=0; $c<=$#group_indices; $c++)
		{
		my @indices = @{$group_indices[$c]};
		my $groupname = shift(@indices);
		my @group_data = @data[@indices];

		my $z_expr_mean=stats::mean(\@group_data);
		my ($x_samples, $x_perc_samples)=(stats::num_above_cutoff(\@group_data, $y_readcount_cutoff));

		# If the absolute number (not percentage) of samples should be used for read count criteria
		if($x_perc_samples_cutoff eq "n")	# If we use '$x_samples_cutoff ne "n"' here we might get some sort of 'argument $x_samples_cutoff not numerical in comparison' error later, because we may unintentionally be forcing perl to interpret $x_samples_cutoff as being a string rather than numeric (?)  
			{
			# If the additional mean read count criterion should NOT be used
			if($z_expr_mean_cutoff eq "n")
				{
				if($x_samples>=$x_samples_cutoff)	{	$pass_abund_criteria++;	}
				}
		
			# Else, if the additional mean read count criterion SHOULD be used 
			else
				{
				if(($x_samples>=$x_samples_cutoff) and ($z_expr_mean>=$z_expr_mean_cutoff))     {       $pass_abund_criteria++; }
				}
			}


		# Else, if the percentage (not absolute number) of samples should be used for abundance criteria
                if($x_samples_cutoff eq "n")
                        {
                        # If the additional mean read count criterion should NOT be used
                        if($z_expr_mean_cutoff eq "n")
                                {
                                if($x_perc_samples>=$x_perc_samples_cutoff)       {       $pass_abund_criteria++; }
                                }

                        # Else, if the additional mean read count criterion SHOULD be used
                        else
                                {
                                if(($x_perc_samples>=$x_perc_samples_cutoff) and ($z_expr_mean>=$z_expr_mean_cutoff))     {       $pass_abund_criteria++; }
                                }
                        }

		} # Loop over groups ends here

	if(($pass_abund_criteria>=$w_groups_cutoff) and ($pass_length eq "yes"))	{	print($out "$line\n");	}
	} # Loop over lines in infile ends here

close($in);
close($out);
print("$proposed_outname");

#close($log);
#close($wlog);
exit;
# end processing

########################################## Define local functions ##########################################

# end functions
