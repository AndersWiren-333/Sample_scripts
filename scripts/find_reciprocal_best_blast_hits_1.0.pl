# DESCRIPTION: This script finds reciprocal best blast hits (RBBHs) by comparing two blast output files (in BLAST outformat 6) representing sets
# of fasta sequences that have been blasted against each other. The goodness of a blast hit is based on bitscore (higher is better), and
# the results are written to a csv file with one line for each RBBH and two columns (sequence id in set 1, sequence id in set 2).
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

# Declare local functions (if any)
sub get_best_hit;

# Define usage (error message to display if the user enters the wrong arguments at the command line)
my $usage = "Usage error for script ${0}. Correct usage: 'perl $0 \$infile1 \$infile2 \$evalue_cutoff \$outname'\n\nwhere".
"\t\$infile1 is the first input file, in BLAST output format 6\n".
"\t\$infile2 is the second input file, in BLAST output format 6\n".
"\t\$evalue_cutoff is the largest e-value a blast hit may have for inclusion\n".
"\t\$outname is the desired name of the outfile that will be produced\n\n";
my @pars = @ARGV or die $usage;
foreach my $el (@pars)  {       $el = text::trim($el);  }
my $infile1 = shift(@pars) or die $usage;
my $infile2 = shift(@pars) or die $usage;
my $evalue_cutoff = shift(@pars) or die $usage;
my $outname = shift(@pars) or die $usage;

# Read infile1 into a matrix
my @matrix1=fileTools::read_table($infile1, "tsv");

# Create best blast hit matrix for @matrix1
my @best_hit_matrix1 = get_best_hit(\@matrix1);

# Free up some memory
@matrix1=();

# Read infile1 into a matrix
my @matrix2=fileTools::read_table($infile2, "tsv");

# Create best blast hit matrix for @matrix1
my @best_hit_matrix2 = get_best_hit(\@matrix2);

# Free up some memory
@matrix2=();


# Compare the two best_hit_matrices to each other and add those records that have each other as best hits to a new matrix
my @reciprocal_matrix=();

HITS1: for(my $c=0; $c<=$#best_hit_matrix1; $c++)
	{
	my $query_1 = $best_hit_matrix1[$c][0];
	my $target_1 = $best_hit_matrix1[$c][1];

	HITS2: for(my $d=0; $d<=$#best_hit_matrix2; $d++)
		{
		my $query_2 = $best_hit_matrix2[$d][0];
		my $target_2 = $best_hit_matrix2[$d][1];

		# Check if the current record in @best_hit_matrix1 matches the current record in @best_hit_matrix2
		if(($query_1 eq $target_2) and ($query_2 eq $target_1))
			{
			push(@reciprocal_matrix, [($query_1, $query_2)]);
			next HITS1;
			}
		}
	}

fileTools::write_table(\@reciprocal_matrix, "csv", $outname, "lin");


#close($log);
#close($wlog);
exit;

# end processing

########################################## Define local functions ##########################################

sub get_best_hit
	{
        my @pars = @_ or die $usage;
        foreach my $el (@pars)  {       $el = text::trim($el);  }
        my $matrixref = shift @pars or die $usage;
	my @matrix = @{$matrixref};

	# Remove records that have an e-value larger than the threshold from the matrix
	# (this actually works, since Perl accepts numbers of the format "2.63e-09" (used in BLAST output files) 
	my @e_matrix=();
	for(my $cc=0; $cc<=$#matrix; $cc++)
		{
		my @arr = @{$matrix[$cc]};
		if($matrix[$cc][10] <= $evalue_cutoff)	{	push(@e_matrix, [(@arr)]);	}
		}

	# Free up some memory
	@matrix=();

	# Sort @e_matrix on query_id (ascending), bit score (descending = larger value better) and e-value (ascending = smaller value better)
	@e_matrix = sort { $a->[0] cmp $b->[0] || $b->[11] <=> $a->[11] || $a->[10] <=> $b->[10] } @e_matrix;

	# Loop over lines in matrix and for each query-id, pick the blast hit with the highest bit score (and the lowest e-value, if there are tied values)
	my @best_hit_matrix=();
	my $comp_qid=333;	# Setting a numerical initial value will give as an error later if we fail to update the value properly (to an alphabetic value), otherwise we may never have noticed if there was an error

	for(my $dd=0; $dd<=$#e_matrix; $dd++)
		{
		my @arr = @{$e_matrix[$dd]};

		# If this is the first line
		if($dd==0)
			{
			$comp_qid = $arr[0];
			push(@best_hit_matrix, [@arr]);	
			}

		# If it is any other line
		else
			{
			# If we are still on the same $query-id, skip this record
			if($arr[0] eq $comp_qid)	{	next;	}

			# If we are on a new query-id
			else
				{
				$comp_qid = $arr[0];
				push(@best_hit_matrix, [@arr]);
				}
			}
		} # Loop over lines in @e_matrix ends
	
	# Free up some memory
	@e_matrix=();

	# Return the best hit matrix
	return(@best_hit_matrix);
	}

# end functions
