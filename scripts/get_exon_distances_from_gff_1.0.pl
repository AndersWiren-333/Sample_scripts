# DESCRIPTION: This script calculates distances between exons in annotated genes, using positional information from a genome annotation file (in gff format). It estimates statistics such
# as mean, median and standard deviation for exon distance, both starnd specifically and overall.
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

# Define usage (error message to display if the user enters the wrong arguments at the command line)
my $usage = "Usage error for script ${0}. Correct usage: 'perl $0 \$infile'\n\nwhere".
"\t\$infile is the annotation file (in gff format) that should be used as input\n".
"\t\$ex_neg is an indicator of whether to exclude genes (RNAs) that have overlapping exons (which produce negative intron lengths. Options are 'y' and 'n')\n\n";
my @pars = @ARGV or die $usage;
foreach my $el (@pars)  {       $el = text::trim($el);  }
my $infile = shift(@pars) or die $usage;
my $ex_neg = shift(@pars) or die $usage;

# Read the matrix
my @matrix=fileTools::read_table($infile, "tsv");

# Loop over matrix lines
my @exons_pos=();
my @exons_neg=();

# Get info about exons on positive and negative strands respectively
LINES: for(my $c=0; $c<=$#matrix; $c++)
	{
	if($matrix[$c][2] eq "exon")
		{
		my $rna=333;
		my $start2="Faster_Tinne";
		my $stop2="Lill_Klas";
		if($matrix[$c][8] =~ /(Parent=)(rna\d+)/)	{	$rna = $2;	}
		else	{	print("Exon has no parent RNA\n"); next LINES;	}
		my $start = $matrix[$c][3];
		my $stop = $matrix[$c][4];
		my $strand = $matrix[$c][6];

		if($start>$stop)
			{
			$start2 = $stop;
			$stop2 = $start;
			$start = $start2;
			$stop = $stop2;
			}

		if($strand eq "+")	{	push(@exons_pos, [$rna, $start, $stop]);	}
		if($strand eq "-")      {       push(@exons_neg, [$rna, $start, $stop]);        }
		}
	}

# Sort exon info matrices on RNA-id, start and stop
@exons_pos = sort { $a->[0] cmp $b->[0] || $a->[1] <=> $b->[1] || $a->[2] <=> $b->[2] } @exons_pos;
@exons_neg = sort { $a->[0] cmp $b->[0] || $a->[1] <=> $b->[1] || $a->[2] <=> $b->[2] } @exons_neg;

my @rnas_pos=();
my @templine=();
my $comp_rna=333;

# Loop over @exons_pos
for(my $d=0; $d<=$#exons_pos; $d++)
	{
	my $rna = $exons_pos[$d][0];
	my $start = $exons_pos[$d][1];
	my $stop = $exons_pos[$d][2];
	my $exon_id = $start."_".$stop;

	# If this is the first line of the matrix
	if($d==0)	{	$templine[0] = $rna; $comp_rna = $rna;	}

	# If it is any other line
	else
		{
		# If we are still on the same RNA
		if($rna eq $comp_rna)
			{
			push(@templine, $exon_id);
			
			# If this is the last line
			if($d==$#exons_pos)	{	push(@rnas_pos, [@templine]);	}
			}

		# If we are on a new RNA
		else
			{
			push(@rnas_pos, [@templine]);
			@templine=();
			$templine[0] = $rna;
			$comp_rna = $rna;
			push(@templine, $exon_id);

                        # If this is the last line
                        if($d==$#exons_pos)     {       push(@rnas_pos, [@templine]);   }
			}
		}
	}

my @rnas_neg=();
@templine=();
$comp_rna=333;

# Loop over @exons_neg
for(my $d=0; $d<=$#exons_neg; $d++)
        {
        my $rna = $exons_neg[$d][0];
        my $start = $exons_neg[$d][1];
        my $stop = $exons_neg[$d][2];
        my $exon_id = $start."_".$stop;

        if($d==0)       {       $templine[0] = $rna; $comp_rna = $rna;  }
        else
                {
                # If we are still on the same RNA
                if($rna eq $comp_rna)
                        {
                        push(@templine, $exon_id);

                        # If this is the last line
                        if($d==$#exons_neg)     {       push(@rnas_neg, [@templine]);   }
                        }

                # If we are on a new RNA
                else
                        {
                        push(@rnas_neg, [@templine]);
                        @templine=();
                        $templine[0] = $rna;
                        $comp_rna = $rna;
                        push(@templine, $exon_id);


                        # If this is the last line
                        if($d==$#exons_neg)     {       push(@rnas_neg, [@templine]);   }
                        }
                }
	}


my @distances_pos=();
my $num_neg_pos=0;

my @distances_neg=();
my $num_neg_neg=0;

# Loop over @rnas_pos
RNAS_POS: for(my $e=0; $e<=$#rnas_pos; $e++)
	{
	my @arr = @{$rnas_pos[$e]};
	if((scalar(@arr)) < 3)	{	next RNAS_POS;	}	# If there is only one exon in the RNA, skip to next RNA
	else
		{
		my $comp_stop="Gundabads_tinnar";

		# Loop over exons in RNA
		for(my $f=1; $f<=$#arr; $f++)
			{
			my ($start, $stop) = split("_", $arr[$f]);
			if($f==1)	{	$comp_stop=$stop;	}
			else
				{
				my $dist = $start-$comp_stop+1;
				if($dist<0)
					{
					print("Negative pos\n");
					$num_neg_pos++;
					if($ex_neg eq "y")	{	next RNAS_POS;	}
					}
				push(@distances_pos, $dist);
				$comp_stop=$stop;
				}
			}
		}
	}


# Loop over @rnas_neg
RNAS_NEG: for(my $g=0; $g<=$#rnas_neg; $g++)
        {
        my @arr = @{$rnas_neg[$g]};
        if((scalar(@arr)) < 3)  {       next RNAS_NEG;   }
        else
                {
                my $comp_stop="Gundabads_tinnar";

                # Loop over exons in RNA
                for(my $h=1; $h<=$#arr; $h++)
                        {
                        my ($start, $stop) = split("_", $arr[$h]);
                        if($h==1)       {       $comp_stop=$stop;       }
                        else
                                {
                                my $dist = $start-$comp_stop+1;
                                if($dist<0)
                                        {
					print("Negative neg\n");
                                        $num_neg_neg++;
                                        if($ex_neg eq "y")      {       next RNAS_NEG;  }
                                        }
                                push(@distances_neg, $dist);
                                $comp_stop=$stop;
                                }
                        }
                }
        }


my $num_neg_all = $num_neg_pos+$num_neg_neg;

# Collect statistics
my @distances_all = (@distances_pos, @distances_neg);

my $n_all=scalar(@distances_all);
my $mean_all=stats::mean(\@distances_all);
my $median_all=stats::median(\@distances_all);
my $percentile_5_all=stats::percentile(\@distances_all, 5);
my $percentile_95_all=stats::percentile(\@distances_all, 95);
my $range_all=stats::range_string(\@distances_all);
my $stdev_all=stats::stdev(\@distances_all);
my $two_stdev_all = 2*$stdev_all;

my $n_pos=scalar(@distances_pos);
my $mean_pos=stats::mean(\@distances_pos);
my $median_pos=stats::median(\@distances_pos);
my $percentile_5_pos=stats::percentile(\@distances_pos, 5);
my $percentile_95_pos=stats::percentile(\@distances_pos, 95);
my $range_pos=stats::range_string(\@distances_pos);
my $stdev_pos=stats::stdev(\@distances_pos);
my $two_stdev_pos = 2*$stdev_pos;

my $n_neg=scalar(@distances_neg);
my $mean_neg=stats::mean(\@distances_neg);
my $median_neg=stats::median(\@distances_neg);
my $percentile_5_neg=stats::percentile(\@distances_neg, 5);
my $percentile_95_neg=stats::percentile(\@distances_neg, 95);
my $range_neg=stats::range_string(\@distances_neg);
my $stdev_neg=stats::stdev(\@distances_neg);
my $two_stdev_neg = 2*$stdev_pos;


# Print statistics to an outfile
open(my $out, ">", "exon_distances.txt") or die "Couldn't create outfile\n";

print($out "Intron length both strands (n = $n_all):\n");
if($ex_neg eq "y")      {       print($out "$num_neg_all negative distances (excluded)\n");     }
else    {       print($out "$num_neg_all negative distances (included)\n");     }
print($out "Mean: $mean_all\n");
print($out "Median: $median_all\n");
print($out "5th percentile: $percentile_5_all\n");
print($out "95th percentile: $percentile_95_all\n");
print($out "Range: $range_all\n");
print($out "Standard deviation: $stdev_all\n");
print($out "2*Standard deviation: $two_stdev_all\n");
print($out "\n\n");

print($out "Intron length positive strand (n = $n_pos):\n");
if($ex_neg eq "y")      {       print($out "$num_neg_pos negative distances (excluded)\n");     }
else    {       print($out "$num_neg_pos negative distances (included)\n");     }
print($out "Mean: $mean_pos\n");
print($out "Median: $median_pos\n");
print($out "5th percentile: $percentile_5_pos\n");
print($out "95th percentile: $percentile_95_pos\n");
print($out "Range: $range_pos\n");
print($out "Standard deviation: $stdev_pos\n");
print($out "2*Standard deviation: $two_stdev_pos\n");
print($out "\n\n");

print($out "Intron length negative strand (n = $n_neg):\n");
if($ex_neg eq "y")      {       print($out "$num_neg_neg negative distances (excluded)\n");     }
else    {       print($out "$num_neg_neg negative distances (included)\n");     }
print($out "Mean: $mean_neg\n");
print($out "Median: $median_neg\n");
print($out "5th percentile: $percentile_5_neg\n");
print($out "95th percentile: $percentile_95_neg\n");
print($out "Range: $range_neg\n");
print($out "Standard deviation: $stdev_neg\n");
print($out "2*Standard deviation: $two_stdev_neg\n");
print($out "\n\n");

#print("Negatives: pos: $num_neg_pos\tneg: $num_neg_neg\n");

close($out);
#close($wlog);
exit;

# end processing

########################################## Define local functions ##########################################

# end functions
