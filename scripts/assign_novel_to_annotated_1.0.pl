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

# Declare local functions (if any)

# Define usage (error message to display if the user enters the wrong arguments at the command line)
my $usage = "Usage error for script ${0}. Correct usage: 'perl $0 \$genelist \$exon_dist_cut \$OFC_cut \$mean_cut \$outname'\n\nwhere".
"\t\$genelist is the file that should be used as input\n".
"\t\$exon_dist_cut\tIf features are closer to each other than this distance (in bp), we consider that they could be exons of the same gene (but it is not the only criterion\n".
"\t\$OFC_cut\tIf the ratio between the offset fold changes for two adjacent features is within this number from 1, we consider that they could be exons of the same gene\n".
"\t\$mean_cut\tIf the ratio between mean gene expression across all samples for two adjacent features is within this number from 1, we consider that they could be exons of the same gene\n".
"\t\$outname is the desired name of the outfile\n\n";

my @pars = @ARGV or die $usage;
foreach my $el (@pars)  {       $el = text::trim($el);  }
my $genelist = shift(@pars) or die $usage;
my $exon_dist_cut = shift(@pars) or die $usage;
my $OFC_cut = shift(@pars) or die $usage;
my $mean_cut = shift(@pars) or die $usage;
my $outname = shift(@pars) or die $usage;

my $ofc_lower = 1-$OFC_cut;
my $ofc_upper = 1+$OFC_cut;
my $mean_lower = 1-$mean_cut;
my $mean_upper = 1+$mean_cut;

# Read genelist
my @matrix=fileTools::read_table($genelist, "csv");

# Extract the header
my $headerref = shift(@matrix);
my @header = @{$headerref};

# a) Find out which column numbers correspond to the sorting keys we will use later
(my $contig_ind) = grep { $header[$_] eq "seqid" } 0..$#header;
(my $start_ind) = grep { $header[$_] eq "start" } 0..$#header;
(my $stop_ind) = grep { $header[$_] eq "end" } 0..$#header;
(my $type_ind) = grep { $header[$_] eq "type" } 0..$#header;
(my $ofc_ind) = grep { $header[$_] =~ /_OFC/ } 0..$#header;
(my $sign_ind) = grep { $header[$_] =~ /_sign/ } 0..$#header;
(my $length_ind) = grep { $header[$_] eq "length" } 0..$#header;

my @data_inds = (1..($length_ind-1));	# Which columns hold the expression data?

# Sort the matrix on contig, start and stop
@matrix = sort { $a->[$contig_ind] cmp $b->[$contig_ind] || $a->[$start_ind] <=> $b->[$start_ind] || $a->[$stop_ind] <=> $b->[$stop_ind] } @matrix;

# Open an outfile to write the results to
open(my $out, ">", $outname) or die "Script $0 couldn't create outfile $outname\n";

# Loop over features in the matrix and check which ones are close enough and similarly expressed enough to belong together
my $parent_gene = 1;
my @new_header = @header;
unshift(@new_header, "Parent_gene");
my $outheader = join(",", @new_header);
print($out "$outheader\n");

my $previous_mean="Graevskopan_Amanda";

for(my $c=0; $c<=$#matrix; $c++)
	{
	my $previous = $c-1;
	my @arr = @{$matrix[$c]};

	my @expression_data = @arr[@data_inds];
	my $mean_expression=stats::mean(\@expression_data);
	
	my $dist = ($arr[$start_ind])-($matrix[$previous][$stop_ind])-1;

	my $ofc_current = $arr[$ofc_ind];
	my $ofc_previous = $matrix[$previous][$ofc_ind];
	if($ofc_current == 0)	{	$ofc_current = 0.01;	}	# We need to be able to for a ratio between the current and previous OFC value, even if one or both are 0. Therefore we add this small offset to any 0 values.
	if($ofc_previous == 0)   {       $ofc_previous = 0.01;    }
	my $ofc_ratio = $ofc_current/$ofc_previous;

	# If this is the first line of the matrix
	if($c==0)
		{
		unshift(@arr, $parent_gene);
		my $outline = join(",", @arr);
		print($out "$outline\n");
		$previous_mean = $mean_expression;
		next;
		}

	# If it is any other line
	else
		{
		my $mean_ratio = $mean_expression/$previous_mean;

		# If this feature is on the same contig as the previous one, is close enough and has a similar expression to it
		if(($arr[$contig_ind] eq $matrix[$previous][$contig_ind]) and ($dist <= $exon_dist_cut) and ($ofc_ratio >= $ofc_lower) and ($ofc_ratio <= $ofc_upper) and ($mean_ratio >= $mean_lower) and ($mean_ratio <= $mean_upper))
			{
			unshift(@arr, $parent_gene);
			my $outline = join(",", @arr);
			print($out "$outline\n");
			}

		# If it isn't or doesn't
		else
			{
			$parent_gene++;
                        unshift(@arr, $parent_gene);
                        my $outline = join(",", @arr);
                        print($out "$outline\n");
			}
		$previous_mean = $mean_expression;
		}
	
 	} # Loop over features in matrix ends

close($out);
#close($log);
#close($wlog);
exit;

# end processing

########################################## Define local functions ##########################################

# end functions
