# DESCRIPTION: No description yet
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

# Define usage (error message to display if the user enters the wrong arguments at the command line)
my $usage = "Syntax error for script ${0}. Correct usage: 'perl $0 \$genelist.csv \$outname'\n\nwhere".
"\t\$genelist is a genelist (in csv format) that should be used as input\n".
"\t\$outname is the desired name of the outfile that will be produced\n\n";

my $infile = shift or die "$usage\n";
my $outname = shift or die "$usage\n";
my $genome = "/scratch2/Anders/01_Projects/00_B_terrestris_genome_files/GCF_000214255.1_Bter_1.0_genomic.fna";
#if(!open(OUT, ">>$outname"))	{	die "$scriptname couldn't create outfile $outname\n";	}

# Read genelist and pick columns
my ($matrix_ref, $header_ref, $col_indices_ref)=misc::file_to_matrix($infile, "comma", "y", "Feature_ID", "type");
my @chosen_ids=();
my @genes = @{$matrix_ref};

# Add top DE genes to @chosen_ids, and strands to @chosen_strands
for(my $c=0; $c<=$#genes; $c++)
	{
	if($genes[$c][1] =~ m/novel/)
		{
		push(@chosen_ids, $genes[$c][0]);
		}
	}

# Get fasta sequences for chosen ids
my @seqs=fastaTools::seqs_from_reference(\@chosen_ids, $genome, $outname);

## Loop over sequences
#for(my $d=0; $d<=$#seqs; $d++)
#	{
#	my $id=$seqs[$d][0];
#	my $seq=$seqs[$d][1];
#	print(OUT ">$id\n$seq\n");
#	}

#close(OUT);
#close(ELOG);
exit;

# end processing

########################################## Define local functions ##########################################

# end functions
