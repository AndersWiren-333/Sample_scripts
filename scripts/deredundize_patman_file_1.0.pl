# DESCRIPTION: This script checks if individual sequence reads align to the genome in several places or not. If they do, they are removed. The infile is a
# tab separated output file from Patman, a sequence aligner.
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
my $usage="Usage: perl check_multiple_matches.pl infile.pat\n";
my $infile=shift or die $usage;
if(!open(IN, $infile))	{	die "Couldn't open $infile\n";	}
if(!open(OUT, ">>nr_$infile"))  {       die "Couldn't create outfile\n";  }
if(!open(STATS, ">>stats_check_multiple.csv"))  {       die "Couldn't create statistics_outfile\n";  }

# Format of line in the patman outfile:
#
#       0       Chromosome/Sequence     NC_015762.1 Bombus terrestris linkage group LG B01, Bter_1.0, whole genome shotgun sequence
#       1       Sequence-abundance      TGAAACAGAAAAAATATCCCATTTAACGAACGCATTTACCGAATCTAGGG-1
#       2       Start                   496
#       3       End                     545
#       4       Strand                  + or -
#       5       No. of mismatches       1

 
# Read patman file into matrix

my @pat_matrix=();

while(<IN>)
	{
	my $line = text::trim($_);
	my @arr = split("\t", $line);
	for(my $i=0; $i<=$#arr; $i++)	{	$arr[$i] = text::trim($arr[$i]);	}
	push(@pat_matrix, \@arr);
	}

my @spat_matrix = sort { $a->[1] cmp $b->[1] } @pat_matrix;
@pat_matrix = @spat_matrix;

# Loop over lines in @pat_matrix

my @comp_arr = @{$pat_matrix[0]};

#foreach my $el (@comp_arr)   {       print("$el\t"); }
#print("\n");


my $comp_seq=$comp_arr[1];
my $hits=1;
my $comp_hits=1;
my $multiple_hit_reads=0;
my $single_hit_reads=0;

for(my $i=1; $i<=$#pat_matrix; $i++)
	{
	$hits++;
	my @arr=@{$pat_matrix[$i]};
	my $seq=$arr[1];

#	print("Comp_seq:\t$comp_seq\nSeq:\t\t$seq\n");

	# If we are still on the same sequence...
	if($seq eq $comp_seq)
		{
		$comp_hits++;	# Increase the hit score for the sequence
		
		# If this is the last line...
		if($i==$#pat_matrix)
			{
			$multiple_hit_reads++;
			}
		}

	else
		{
		# If hit score is 1...
		if($comp_hits == 1)
			{
			# Print previous sequence (= @comp_arr) to the outfile
			$single_hit_reads++;
			my $outline=join("\t", @comp_arr);
			print(OUT "$outline\n");
			
			# Set new initials
			@comp_arr=@arr;
			$comp_seq=$seq;
			$comp_hits=1;
	
			# If this is the last line...
			if($i==$#pat_matrix)
				{
				$single_hit_reads++;
	                        my $outline=join("\t", @comp_arr);
                        	print(OUT "$outline\n");
				}

			}
	
		# If not, the comp_seq is a multiply aligning sequence read and should not be written to the outfile 
		else
			{
			$multiple_hit_reads++;

			# Set new initials
                        @comp_arr=@arr;
                        $comp_seq=$seq;
                        $comp_hits=1;

                        # If this is the last line...
                        if($i==$#pat_matrix)
                                {
				$single_hit_reads++;
                                my $outline=join("\t", @comp_arr);
                                print(OUT "$outline\n");
                                }
			}
		}
	}

my $total=$multiple_hit_reads+$single_hit_reads;
print(STATS "$infile,$hits,$total,$single_hit_reads,$multiple_hit_reads\n");

close(IN);
close(OUT);
close(STATS);
exit;

# end processing

########################################## Define local functions ##########################################

# end functions
