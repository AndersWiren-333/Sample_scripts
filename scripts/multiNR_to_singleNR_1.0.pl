# DESCRIPTION: This script concatenates a set of non-redundant fasta files, one by one, and makes the concatenated file non-redundant.
# A single input file is a list of unique sequence reads with the abundance of each read in that sample, and each read
# only occurs once in each file. However, the same read may occur in two different files. If it does, after those two files
# have been concatenated, the script merges the list entry for that read from two into one, summing the abundances of the read.
# The script starts by merging the first two input files in this way, and then merges the new (probably larger) file with the
# next input file in a list of all files that should be merged. By doing the merge in this stepwise fashion, the process will
# require less memory than if all files were first concatenated and then had their read entires merged. 
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

# Declare variables and filehandles
my $usage="perl multiNR_to_singleNR.pl filenames_file\n";		# filenames_file is a textfile with the names of the files that should be concatenated (one filename on each row)
my $namefile=shift or die $usage;
if(!open(NAM, $namefile))	{	die "Couldn't open namefile\n";	}
if(!open(LOG, ">>logfile_redund.txt"))	{	die "Couldn't create logfile\n";	}

# Make a list of infiles
my @files=();
while(<NAM>)
	{
	my $filename=$_;
	$filename=text::trim($filename);
	push(@files, $filename);
	}
close(NAM);


# Concatenate and merge infiles one by one
my @grand_matrix=();				# Create matrix to hold concatenated and merged data
my @new_matrix=();
for(my $i=0; $i<scalar(@files); $i++)	# Loop over files
	{
	my $file=$files[$i];
	my $fa = Bio::SeqIO->new('-format'=>'fasta', '-file' => "$file");

	# Read in data and push at matrix
	while(my $seq = $fa -> next_seq())
		{
		my $sequence = $seq->seq();
		my $id = $seq->display_id();

		my @id_parts=split("-", $id);
		my $new_id=$id_parts[0];
		my $abund=$id_parts[1];

		push(@grand_matrix, [$sequence, $abund]);
		}
	print(LOG "Finished adding file $i to matrix\n");

	# If it's file 1, don't merge anything
	if($i==0)	{	print(LOG "File $i finished\n");	}

	# Else, merge
	else
		{
		@grand_matrix = sort { $a->[0] cmp $b->[0] } @grand_matrix;	# sort matrix
		print(LOG "Finished sorting matrix after file $i\n");

		# Merge reads

		# Set initials
		my $current_seq=$grand_matrix[0][0];
		my $abu_count=$grand_matrix[0][1];


		for(my $k=1; $k<scalar(@grand_matrix); $k++)	# Loop over reads...
			{
			if($k != $#grand_matrix)	# As long as it isn't the last read in the matrix...
				{
				if($grand_matrix[$k][0] eq $current_seq)        {       $abu_count=$abu_count+$grand_matrix[$k][1];     }       # if it is a copy of the previous read, add abundance
				else
					{
					push(@new_matrix, [$current_seq, $abu_count]);  # if it is a new read, add previous read (and its abundance) to new matrix
                                        $current_seq=$grand_matrix[$k][0];      # set new initial read
                                        $abu_count=$grand_matrix[$k][1];        # and new initial abundance
					}
				}		
			else			# If you have come to the last read...
				{
				if($grand_matrix[$k][0] eq $current_seq)	{	$abu_count=$abu_count+$grand_matrix[$k][1]; push(@new_matrix, [$current_seq, $abu_count]);}
				else
				        {
                                        push(@new_matrix, [$current_seq, $abu_count]);  # if it is a new read, add previous read (and its abundance) to new matrix
                                        $current_seq=$grand_matrix[$k][0];      # set new initial read
                                        $abu_count=$grand_matrix[$k][1];        # and new initial abundance
					push(@new_matrix, [$current_seq, $abu_count]);	# add the final line
                                        }

				}
			}
		
		@grand_matrix=@new_matrix;
		@new_matrix=();
		print(LOG "Finished merging after file $i\n");
		}
	}

if(!open(OUT, ">>nr_all.fa"))	{	die "Couldn't create outfile";	}
print(LOG "Now writing to outfile...");
for(my $j=0; $j<scalar(@grand_matrix); $j++)
	{
	my $sek=$grand_matrix[$j][0];
	my $abu=$grand_matrix[$j][1];
	my $id_line=">$sek-$abu";
	print(OUT "$id_line\n$sek\n");
	}
print(LOG " Finished!");

close(LOG);
close(OUT);
exit;

# end processing

########################################## Define local functions ##########################################

# end functions
