# DESCRIPTION:
# This script picks a number (given at command line) of sequences at random from a pair of paired-end fastq files
# and prints them to two corresponding fastq outfiles.

#	NB NB!!
#	The files have to be synchronised, i.e. the reads in the "left" (or "R1") file need to be
#	listed in the same order as the reads in the "right" (or "R2") file. This can be achieved by using 
# 	pe_sync_fastq_1.0.pl.
#	NB NB!!

# If the original files contain more than a given number of sequences, the script divides them into parts (of that many
# sequences each) and subsamples the parts instead (otherwise it may be terribly slow). I recommend 3 million sequences per part
# for best performance. The end result is the same, two fastq files with the subsample. The two original files are retained.

# This is how to use the script from the command line: "perl pe_subsample.pl seqs_per_part to_num file_01 file_02"

	# seqs_per_part = number of sequences per part that infiles are (temporarily) split into
	# to_num = number of sequences to pick out
	# file_01 = infile 1 (the "left" or "R1" file)
	# file_02 = infile 2 (the "right" or "R2" file)
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
my $scripts = cwd;
chdir("../modules");
my $modules = cwd;
chdir("../scripts");

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
sub divide_fastq;
sub read_fastq;
sub testprint_matrix_table;
sub testprint_matrix_list;
sub print_matrix;
sub transpose_matrix;
sub deredund_numeric_list;

# Declare variables and filehandles
my $usage = "Correct syntax: perl pe_subsample.pl seqs_per_part to_num file_01 file_02";
my $seqs_per_part = shift or die $usage;
my $sample_to_num = shift or die $usage;
my $file_01 = shift or die $usage;
my $file_02 = shift or die $usage;

# Divide first file into parts
my @parts_and_last_1 = divide_fastq($file_01, $seqs_per_part, "r1");
my $num_parts = $parts_and_last_1[0];
my $total_seqs = $parts_and_last_1[1];
my $last_seqs = $parts_and_last_1[2];

# Divide first file into parts
my @parts_and_last_2 = divide_fastq($file_02, $seqs_per_part, "r2");

# How many reads should be picked from each part? And how many should be picked from the last (smaller) part?
my $prop_whole_part = $seqs_per_part/$total_seqs;
my $prop_last_part = $last_seqs/$total_seqs;
my $sample_to_whole_part = int(${prop_whole_part}*$sample_to_num);
my $sample_to_last_part = int(${prop_last_part}*$sample_to_num);

# Adjust the number of sequences to choose (the number will ususally be somewhat underestimated by the rounding of numbers above)
my $total_chosen = (($num_parts-1)*($sample_to_whole_part))+$sample_to_last_part;
my $extra = (${sample_to_num}-${total_chosen});

# Chose the parts to pick the extra sequences from
my @special_parts=();
for(my $e=1; $e<=$extra; $e++)
	{
	my $part = (int(rand($num_parts)))+1;
	print("I have chosen file $part\n");
	push(@special_parts, $part);
	}


my @special_scores = deredund_numeric_list(\@special_parts);

print("Sample to $sample_to_num\n");
print("Total chosen: $total_chosen\n");
print("We need $extra more sequences\n");


print("\n-----------------\n");
testprint_matrix_table(\@special_scores);


# Loop over (pairs of) parts
for(my $d=1; $d<=$num_parts; $d++)
	{
	my $to=0;
	my $from=0;
	my $extra=0;

	# Check if this part is among those chosen to pick extra sequences from
	# Loop over @special_scores
	for(my $i=0; $i<=$#special_scores; $i++)
		{
		if($d==$special_scores[$i][0])	{	$extra = $special_scores[$i][1];	}
		}

	# Set different sample to and sample from numbers depending on whether this is a whole or last part
	if($d != $num_parts)	{	$to = $sample_to_whole_part+$extra; $from = $seqs_per_part;	}
	elsif($d == $num_parts)    {       $to = $sample_to_last_part+$extra; $from = $last_seqs;   }

	# Create the correct filenames to open the current part files
	my $first = "r1_".$d;
	my $second = "r2_".$d;

	# Read files into matrices
	my @one = read_fastq($first);
	my @two = read_fastq($second);

	# Declare matrices to hold chosen reads
	my @chosen_seqs_1=();
	my @chosen_seqs_2=();

	print("Files $file_01 and $file_02, part $d. Choosing $to sequences out of ${from}...\n");

	# Pick $to number of sequences at random (loop $to times)
	for(my $c=1; $c<=$to; $c++)
		{
		my $num = int(rand($from));		# The number drawn can be 0, but not $from (but that's alright, because the last index of the fastq matrices are $from-num -1)
							# So we are choosing a matrix index.

		# Add the chosen sequence to new outfile matrix and remove it from the old matrix
		push(@chosen_seqs_1, $one[$num]);
		splice(@one, $num, 1);

		# Do the same for the other fastq matrix
        	push(@chosen_seqs_2, $two[$num]);
        	splice(@two, $num, 1);

		# Decrease the number of sequences to choose randomly from in next cycle by one (since there is now one less sequence to choose from)
		$from--;
		}

	# Print chosen sequences to outfiles
	print_matrix(\@chosen_seqs_1, "sub_$file_01");
	print_matrix(\@chosen_seqs_2, "sub_$file_02");

	# Delete temporary files
	unlink($first, $second);
	}

exit;

# end processing

########################################## Define local functions ##########################################

sub divide_fastq
	{
	my $file = text::trim(shift);
	my $limit = text::trim(shift);
	my $end = text::trim(shift);


	# Open the file and write sequences to outfiles
	if(!open(IN, $file))	{	die "Couldn't open infile $file for division into parts\n";	}
	my $lineno = 1;
	my $part = 1;
	my $limit_base = 4*$limit;
	$limit = $limit_base;
	my $total_lines=0;
	if(!open(OUT, ">>${end}_${part}"))	{	die "Couldn't create intermediate file ${end}_${part} for infile $file\n";	}
	#print("Created file ${end}_${part}\n");

	while(<IN>)
		{
		my $line = text::trim($_);
		$total_lines++;
		if($lineno<$limit)	{	print(OUT "$line\n");	}
		elsif($lineno==$limit)
			{
			print(OUT "$line\n");
			close(OUT);
			$limit=$limit+$limit_base;
			$part++;
			if(!open(OUT, ">>${end}_${part}"))      {       die "Couldn't create intermediate file ${end}_${part} for infile $file\n";    }
			#print("Created file ${end}_${part}\n");
			}
		$lineno++;
		}
	close(OUT);
	close(IN);
	my $total_seqs = $total_lines/4;

	# Open the last part and count the number of sequences in it (it will most likely have less sequences than the others)
	if(!open(LAST, "${end}_${part}"))	{	die "Couldn't open last part of infile $file for counting sequences\n";	}	
	my $last_lines=0;
	while(<LAST>)	{	$last_lines++;	}
	close(LAST);

	if($last_lines==0)
		{
		#print("\tNo lines in last part. Must be an even number of sequences :-)\n");
		$last_lines = $limit_base;
		unlink("${end}_${part}");
		#print("Deleted file ${end}_${part}\n");
		$part--;
		}

	#print("Last lines is $last_lines\n");
	my $number_of_parts = $part;
	my $last_seqs = $last_lines/4;
	my @parts_and_last = ($number_of_parts, $total_seqs, $last_seqs);

	#print("Last part is $part\n");
	return(@parts_and_last);
	}

sub read_fastq
        {
        my $file = shift;
        if(!open(FIL, $file))   {       die "Couldn't open $file_01\n"; }
        my @matrix=();
        my @temp=();
        my $counter=0;
        while(<FIL>)
                {
                $counter++;
                my $line = text::trim($_);
                push(@temp, $line);
                if($counter==4)                         # If this is the last line of a sequence record...
                        {
                        push(@matrix, [@temp]);         # Add the four lines of the record as an element (a line of) to @one matrix
                        @temp=();                                       # Free the @temp array
                        $counter=0;                             # Reset the counter
                        }
                }
        return(@matrix);
        }

sub testprint_matrix_table
        {
        my $matrix_ref = shift;
        my @matrix = @{$matrix_ref};
        for(my $c=0; $c<=$#matrix; $c++)
                {
                my @arr=@{$matrix[$c]};
                foreach my $el (@arr)   {       print("$el\t"); }
		print("\n");
                }
        }

sub testprint_matrix_list
        {
        my $matrix_ref = shift;
        my @matrix = @{$matrix_ref};
        for(my $c=0; $c<=$#matrix; $c++)
                {
                my @arr=@{$matrix[$c]};
                foreach my $el (@arr)   {       print("$el\n"); }
                }
        }


sub print_matrix
        {
        my $matrix_ref = shift;
        my $outname = shift;
        my @matrix = @{$matrix_ref};
        if(!open(OUT, ">>$outname"))    {       die "Couldn't create outfile $outname\n";       }
        for(my $c=0; $c<=$#matrix; $c++)
                {
                my @arr=@{$matrix[$c]};
                foreach my $el (@arr)   {       print(OUT "$el\n");     }
                }
        close(OUT);
        }

sub transpose_matrix
        {
        my $matrix_ref = shift;
        my @matrix = @{$matrix_ref};
        my $num_rows=scalar(@matrix);
        my @first_row=@{$matrix[0]};
        my $num_cols=scalar(@first_row);

        my @t_matrix=();

        # Loop over columns in original matrix
        for(my $cc=0; $cc<$num_cols; $cc++)
                {
                # Loop over rows in riginal matrix
                for(my $dd=0; $dd<$num_rows; $dd++)
                        {
                        $t_matrix[$cc][$dd] = $matrix[$dd][$cc];
                        }
                }

        return(@t_matrix);
        }

sub deredund_numeric_list
	{
	my $arref = shift;
	my @arr = @{$arref};
	@arr = sort { $a <=> $b } @arr;

	my @out_matrix=();
	my $ref = $arr[0];
	my $ref_score = 1;

	for(my $cc=1; $cc<=$#arr; $cc++)
		{
		my $value = $arr[$cc];
		if($cc==$#arr)
			{
			if($value==$ref)
				{
				$ref_score++;
				push(@out_matrix, [($ref, $ref_score)]);
				}
			else
				{
				push(@out_matrix, [($ref, $ref_score)]);
				$ref=$value;
				$ref_score=1;
				push(@out_matrix, [($ref, $ref_score)]);
				}
			}
		
		else
			{
			if($value==$ref)	{	$ref_score++;	}
			else
				{
				push(@out_matrix, [($ref, $ref_score)]);
				$ref=$value;
				$ref_score=1;
				}
			}
		} # Loop over array stops here

	# A line in the new matrix has the format:
	
	#	0	number					5
	#	1	times that number occures in the list	10

	testprint_matrix_table(\@out_matrix);
	return(@out_matrix);
	} # Function ends here

# end functions
