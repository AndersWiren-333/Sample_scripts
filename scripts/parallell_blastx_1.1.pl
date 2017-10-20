# DESCRIPTION: This script reads a set fasta files in parallell. It sends one sequence at the time from a file to the NCBI to be BLASTx:ed against the nr database. The output (blast hits, max 100 per sequence)
# is written to one textfile for each fasta infile.
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

# Set paths to scripts and modules. Setting explicit paths to the scripts and modules in this specific repository (rather than adding paths to @INC, PERLLIB and PATH on
# your local system) avoids the risk of scripts calling the wrong scripts/modules if you have other repositories on your system that happen to have some script- and module names
# in common with this repository.
my $thisfile = (__FILE__);
my $scripts = "";
if($thisfile =~ m/^(.+)\//)	{	$scripts = $1;	}
my $modules = $scripts;
$modules =~ s/scripts/modules/;

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
sub run_blast;

# Declare variables and filehandles
my $usage = "perl parallell_blastx.pl folder_with_infiles number_of_parts_to_split_the_job_into";
my $folder = shift or die $usage;
my $num_threads = shift or die $usage;
$folder="$folder\\";
my @threads=(1..$num_threads);

foreach(@threads)	{	$_ = threads->create(\&run_blast);	}
foreach(@threads)	{	$_->join();	}

print("\n\n\tFinally finished!\n\n");
exit;

# end processing

########################################## Define local functions ##########################################

sub run_blast
	{
	# Define variables
	my $grand_start=time();
	my $id = threads->tid();
	my $infile = ($folder . "${id}.fa");
	my $outfile=($folder . "blast_results_${id}" . ".txt");
	my $temp_in=($folder . "temp_in_${id}.fa");
	my $temp_out=($folder . "temp_out_${id}.txt");
	my $log=($folder . "blast_logfile_${id}.txt");
	my $redo=($folder . "reblast_these_${id}.fa");
	
	# Open filehandles
	if(!open(IN, $infile))	{	die "Couldn't open infile $infile";	}
	if(!open(OUT, ">>$outfile"))	{	die "Couldn't create outfile_${id}.txt\n";	}
	if(!open(LOG, ">>$log"))	{	die "Couldn't create logfile for file_${$id}\n";	}
	if(!open(REDO, ">>$redo"))	{	die "Couldn't create reblast file for file_${$id}\n";	}

	# Loop over lines in the infile
	my $seq_index=1;	
	my $header="";
	my $seq="";
	my $time_so_far = 0;
	my $mean_time = 0;
	
	while(<IN>)
		{
		my $line = trim($_);		# Read in line
		if($line =~ />/)			# If it is the header line of a sequence...
			{
			$header=$line;
			if(!open(TEMP_IN, ">>$temp_in"))	{	die "Couldn't create tempfile\n";	}	# ... open a temporary infile
			print(TEMP_IN "$line\n");			# print the header line to temporary infile (the file we will send to blast later)
			}
		else		# If it is the sequence line of the fasta sequence...
			{
			my $start_time = time();
			$seq=$line;
			print(TEMP_IN "$line");															# print the sequence line to the temporary infile (the file we will send to blast later)
			close(TEMP_IN);																		# close the temporary infile (blast will open it by itself later)
			
			# Put together the blast command
			my $first="\"..\\..\\..\\..\\Program Files\\NCBI\\blast-2.5.0+\\bin\\blastx.exe\" -db nr -query ";
			my $second=" -outfmt \"6 qseqid stitle evalue bitscore score length pident\" -max_target_seqs 100 -evalue 1e-005 -query_gencode 1 -word_size 6 -matrix \"BLOSUM62\" -gapopen 11 -gapextend 1 -comp_based_stats 2 -seg \"yes\" -remote";
			my $cmd=($first . $temp_in . " -out " . $temp_out . $second);

			# Run the command and capture potential error messages
			print("Now blasting file $id sequence $seq_index...\n");
			my $result=`$cmd 2>&1`;
			my $exit_code=$?;
			my $stop_time = time();
			my $total_time = ($stop_time-$start_time)/60;
			$time_so_far = $time_so_far+$total_time;
			$mean_time = ($time_so_far/$seq_index);
			
			if($result ne "")
				{
				print(LOG "$header\t$result\t(exit code $exit_code)\n------------------\n");
				print(REDO "$header\n$seq\n");
				}
			
			# Write results from temporary outfile to permanent outfile
			if(!open(TEMP_OUT, "$temp_out"))	{	die "Couldn't open temporary results file out_temp.txt\n";	}
			my @results=();
			while(<TEMP_OUT>)
				{
				my $line=trim($_);
				push(@results, $line);
				}
				
			for(my $i=0; $i<=$#results; $i++)
				{
				print(OUT "$results[$i]\n");
				}
			close(TEMP_OUT);
			
			# Clean up
			print("Finished file $id sequence $seq_index\n");
			print(LOG "Finished $header (sequence $seq_index) (total time: $total_time min,  mean time now: $mean_time min)\n");
			$seq_index++;
			unlink($temp_in);
			unlink($temp_out);
			}
		}
	my $grand_stop=time();
	my $run_time = ($grand_stop-$grand_start)/3600;
	
	print(LOG "\nRuntime (hours):\t$run_time\n");
	close(IN);
	close(OUT);
	close(LOG);
	close(REDO);
	}

# end functions
