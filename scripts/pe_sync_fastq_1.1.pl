# DESCRIPTION: This script synchronizes paired end fastq files. It takes two infiles, e.g. sample_A_R1.fastq (containing the "left" reads) and sample_A_R2.fastq (containing the "right" reads).
# The output is two new files, sync_sample_A_R1.fastq and sync_sample_A_R2.fastq that contain only reads from the R1 file that have a corresponding read in the R" file and vice versa.
# Unpaired reads (= singletons = reads only found in one of the two infiles) are written to two separate outfiles (singletons_sample_A_R1.fastq and singletons_sample_A_R2.fastq).
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

# Define usage (error message to display if the user enters the wrong arguments at the command line)
my $usage = "Syntax error for script ${scriptname}. Correct usage: 'perl $scriptname \$file_R1.fastq \$file_R2.fastq'\n\nwhere".
"\t\$file_R1.fastq is a fastq file containing the 'left' reads\n".
"\t\$file_R2.fastq is a fastq file containing the 'right' reads\n\n";

my $file_01 = shift or die $usage;
my $file_02 = shift or die $usage;
my $newname_01 = "sync_$file_01";
my $newname_02 = "sync_$file_02";
my $sing_01 = "singletons_$file_01";
my $sing_02 = "singletons_$file_02";

# Read files into matrices
my @one = fastqTools::read_fastq($file_01);
my @two = fastqTools::read_fastq($file_02);

# Format of a line in the fastq matrices:

#	0	Sequence ID		@D1SF08P1:78:C7RTVACXX:4:1101:10000:100055 1:N:0:ATGTCA
#	1	Sequence		TCTTCATTTAAAATTTATCTTAATTCAACATCGAGGTCGCAATCATCTTTTTCAATAAGATCTTTAAAAAGATATTACGCTGTTATCCCTAAGGTAATTTA
#	2	Optional description	+
#	3	Quality scores		CCCFFFFFGHGHHJJJJJJJJJJJJJJJJJJJIJJFGHIJJJJJJJJJJJGIIJIJJJIJJIJIJJEIGIGHHGHEHFFFFDDEDEDDDDDCDC>ACDEEC

# Sort matrices alphabetically on id
@one = sort { $a->[0] cmp $b->[0] } @one;
@two = sort { $a->[0] cmp $b->[0] } @two;

# Declare matrices to hold singletons
my @singletons_01=();
my @singletons_02=();

# Compare matrices
my $counter_01=0;
my $counter_02=0;
my @new_01=();
my @new_02=();

# Declare id variables
my $id_01="";
my $id_02="";

# Loop over matrices in parallell, as long as the comp_id's are defined
MATLINES: while((defined($one[$counter_01][0])) and (defined($two[$counter_02][0])))
	{
	# Set id's (if we have gone to the end of either file, exit loop)
	if($one[$counter_01][0] =~ /(@\w+:\d{2}:\w+:\d:\d+:\d+:\d+) (\d:\w:\d:\w{6})/)	{	$id_01=$1;	}
	if($two[$counter_02][0] =~ /(@\w+:\d{2}:\w+:\d:\d+:\d+:\d+) (\d:\w:\d:\w{6})/)	{	$id_02=$1;	}

	# Compare id's
        if($id_01 eq $id_02)
                {
                push(@new_01, $one[$counter_01]);
                push(@new_02, $two[$counter_02]);
                $counter_01++;
                $counter_02++;
                }
	elsif($id_01 lt $id_02)
		{
		push(@singletons_01, $one[$counter_01]);		# Add that sequence to the singletons list for matrix @one
		$counter_01++;
		}
	elsif($id_01 gt $id_02)
		{
		push(@singletons_02, $two[$counter_02]);                # Add that sequence to the singletons list for matrix @two
		$counter_02++;
		}
	}

my $ett = scalar(@new_01);
my $tva = scalar(@new_02);

print("Sequences in new file 1: $ett\n");
print("Sequences in new file 2: $tva\n");

fastqTools::print_fastq(\@new_01, $newname_01);
fastqTools::print_fastq(\@new_02, $newname_02);
fastqTools::print_fastq(\@singletons_01, $sing_01);
fastqTools::print_fastq(\@singletons_02, $sing_02);

exit;

# end processing

########################################## Define local functions ##########################################

# end functions
