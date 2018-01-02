# DESCRIPTION: Calculates sample complexities when succesively subsampling a set of fasta files to smaller total
# number of sequences.
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
my $usage = "perl run_resampling_check infile_1 infile_2 infile_3 etc...\n";

my @samples=qw(
test_01.fa.gz
test_02.fa.gz

test_03.fa.gz
test_04.fa.gz

test_05.fa.gz
test_06.fa.gz

test_07.fa.gz
test_08.fa.gz);

=pod
my @samples=qw(
rRa_nr_01_B24_1.fa.gz
rRa_nr_02_B24_2.fa.gz
rRa_nr_03_B24_3.fa.gz
rRa_nr_04_B96_1.fa.gz
rRa_nr_05_B96_2.fa.gz
rRa_nr_06_B96_3.fa.gz
rRa_nr_07_B96_Low.fa.gz

rRa_nr_08_BLay_1.fa.gz
rRa_nr_09_BLay_2.fa.gz
rRa_nr_10_BLay_3.fa.gz
rRa_nr_11_O24_1.fa.gz
rRa_nr_12_O24_2.fa.gz
rRa_nr_13_O24_3.fa.gz
rRa_nr_14_O96_1.fa.gz

rRa_nr_15_O96_2.fa.gz
rRa_nr_16_O96_3.fa.gz
rRa_nr_17_OLay_1.fa.gz
rRa_nr_18_OLay_2.fa.gz
rRa_nr_19_OLay_3.fa.gz
rRa_nr_20_FB24_1.fa.gz
rRa_nr_21_FB24_2.fa.gz

rRa_nr_22_FB24_3.fa.gz
rRa_nr_23_FB96_1.fa.gz
rRa_nr_24_FB96_2.fa.gz
rRa_nr_25_FB96_3.fa.gz
rRa_nr_26_FBLay_1.fa.gz
rRa_nr_27_FBLay_2.fa.gz
rRa_nr_28_FBLay_3.fa.gz);
=cut



# Set absolute path to scripts folder
my $scripts="/scratch/Anders/02_test/01_scripts";
my $ref_genome="/scratch/Anders/02_test/01_scripts/GCF_000214255.1_Bter_1.0_genomic.fna";

#my $scripts="/gpfs/home/zjh15phu/01_Scripts/06_Resampling_check";
#my $ref_genome="/gpfs/home/zjh15phu/01_Scripts/06_Resampling_check/GCF_000214255.1_Bter_1.0_genomic.fna";
my $compress="gzip";
my $uncompress="gunzip";

my @numbers=(0.95, 0.90, 0.85, 0.80, 0.75, 0.70, 0.65, 0.60, 0.55, 0.50);
my $stat_header = "File,Total_reads,Unique_reads,Complexity,Copies_per_read,Total_genome_matching_reads_1mm,Perc_genome_matching_reads_1mm,".
"Unique_genome_matching_reads_1mm,Perc_unique_genome_matching_reads,Complexity_genome_matching_reads_1mm,Copies_per_read_genome_matching_reads_1mm";
if(!open(STAT, ">>stats_subsampling_check.csv"))	{	die "Couldn't create statistics outfile\n";	}
print(STAT "$stat_header\n");

# Loop over samples
for(my $c=0; $c<=$#samples; $c++)
	{
	my $nr_fasta_gz = $samples[$c];
	my $nr_fasta = $nr_fasta_gz;
	$nr_fasta =~ s/.gz//;
	system("$uncompress $nr_fasta_gz");
	
	# Subsample to a range of percentages (loop over @numbers)
	for(my $d=0; $d<=$#numbers; $d++)
		{
		my @stats_vector=();			# Set up an array to hold statistics for this particular subsampled file
		my $perc = $numbers[$d];
		my $r_sub_fasta = "r_${perc}_${nr_fasta}";

		my $patfile = "frigg.pat";
		#$patfile =~ 's/fa/pat/';

		# Subsample
		system("perl $scripts/bootstrap.pl $nr_fasta $perc -1 $r_sub_fasta");
		system("perl $scripts/R2NR.pl $r_sub_fasta frigg.fa");

		# Count reads in resulting file
		my $reads = `perl $scripts/count_reads_in_nr_fasta_1.0.pl frigg.fa`;
		$reads = text::trim($reads);
		my @read_counts = split(",", $reads);
		my $total_reads = $read_counts[0];
		my $uni_reads = $read_counts[1];
		my $complexity=($uni_reads/$total_reads);
		my $cop_per_read=1/$complexity;
		push(@stats_vector, $r_sub_fasta, $total_reads, $uni_reads, $complexity, $cop_per_read);

		# Align to genome
		system("patman -D $ref_genome -P frigg.fa -o $patfile -e 1");
		unlink($r_sub_fasta);
		unlink("frigg.fa");

		# Count genome matching reads
		my $read_counts_string = `perl $scripts/count_genome_matching_reads_2.1.pl $patfile`;
		$read_counts_string = text::trim($read_counts_string);
		@read_counts = split(",", $read_counts_string);
		my $gen_total_reads = $read_counts[0];
		my $gen_uni_reads = $read_counts[1];
		my $perc_gen_total_reads = ($gen_total_reads/$total_reads)*100;
		my $perc_gen_uni_reads = ($gen_uni_reads/$uni_reads)*100;
        	my $gen_complexity = ($gen_uni_reads/$gen_total_reads);
		my $gen_cop_per_read = 1/$gen_complexity;
        	push(@stats_vector, $gen_total_reads, $perc_gen_total_reads, $gen_uni_reads, $perc_gen_uni_reads, $gen_complexity, $gen_cop_per_read);
		my $out_stats=join(",", @stats_vector);
		print(STAT "$out_stats\n");
		unlink($patfile);

		} # Loop over percentages ends here

	system("$compress $nr_fasta");	

	} # Loop over samples ends here

close(STAT);
exit;

# end processing

########################################## Define local functions ##########################################

# end functions
