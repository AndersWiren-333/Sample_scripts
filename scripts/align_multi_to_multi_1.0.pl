# DESCRIPTION: This script aligns multiple files with sequence reads (in fasta format) to multiple target files in fasta format
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

# Define usage (error message to display if the user enters the wrong arguments at the command line)
my $usage = "Syntax error for script ${0}. Correct usage: 'perl $0 \$query_names.txt \$target_names_and_mismatches.txt \$fasta_format \$zip'\n\nwhere".
"\t\$query_names.txt is a textfile listing the names (one on each line) of the non-redundant fasta files that should be aligned to something else\n".
"\t\tThe first line needs to be the absolute path to the folder that the files are in, and all other lines are the names of the files\n".
"\t\$target_names_and_mismatches.txt is a textfile listing each fasta file (absolute path if necessary) that the query files should be aligned to, and the number of mismatches that should be used in alignment\n".
"\t\for example:\n".
"\t\t'reference_genome_1.fa,0,1'\n".
"\t\t'collection_of_miRNAs1.fa,0,1,2'\n".
"\t\$fasta_format is the format of the fasta files to be aligned. This can be either 'r', 'nr' or 'mnr'\n\nwhere".
"\t\t'r' is 'redundant' i.e. '>id (linebreak) sequence'\n".
"\t\t'nr' is 'non-redundant' i.e. '>sequence-abundance (linebreak) sequence'\n".
"\t\t'mnr' is 'multi-non-redundant' i.e. '>sequence_abundance1_abundance2_abundance3... (linebreak) sequence'\n".
"\t\$zip is an indicator of whether either input or output files or both should be compressed after the script has finished (options are 'in', 'out' and 'all'. Default is no compression)\n";

# Accept input parameters and open filehandles
my @pars = @ARGV or die $usage;
foreach my $el (@pars)	{	$el = text::trim($el);	}
my $queries = shift(@pars) or die $usage;
my $targets = shift(@pars) or die $usage;
my $fasta_format = shift(@pars) or die $usage;
my $zip = "no";
if(defined($pars[0]))	{	$zip = $pars[0];	}
my $script_basename = (text::parse_path($0))[2];

my $statname = "stats_${script_basename}_${timestamp}.csv";
open(my $stat, ">>", $statname) or die "Script $0 couldn't create statistics outfile\n";

# Read the query names
my @query_files=misc::read_list($queries);
my $query_folder = shift(@query_files);

# Read the target names and their respective mismatch numbers
my @raw_targets=misc::read_list($targets);
my @targets=();

# Create a header for the statistics outfile
my @header1=();
my @header2=();
push(@header1, "Reference ->", "", "", "");
push(@header2, "File", "", "Total_reads", "Unique_reads");

# Loop over targets to create headers
for(my $c=0; $c<=$#raw_targets; $c++)
	{
	my @arr = split(",", $raw_targets[$c]);
	push(@targets, \@arr);
	my $target_raw = $arr[0];
	my @parts = split("/", $target_raw);
	if($parts[$#parts] eq "/")	{	pop(@parts);	}
	my $target_name="$parts[$#parts]";

	# Loop over mismatches
	for(my $z=1; $z<=$#arr; $z++)
		{
		my $mism = $arr[$z];
		my $header1_name = "$target_name"." (${mism} mismatches)";
		push(@header1, "", $header1_name, "", "", "", "", "", "");
		push(@header2, "", "total_matching_reads", "perc_tot_match_reads", "unique_matching_reads", "perc_uni_match_reads", "complexity", "copies_per_read", "multi_matching_reads");	
		}
	}

my $out_header1 = join(",", @header1);
my $out_header2 = join(",", @header2);
print($stat "$out_header1\n$out_header2\n");

# Loop over query files
for(my $d=0; $d<=$#query_files; $d++)
	{
	my @stats_vector=();

	my $query_path = "${query_folder}/"."$query_files[$d]";
	$query_path =~ s/.gz//;
	my $query_file = $query_files[$d];
	
	# Uncompress query file if necessary
	if($query_file =~ m/.gz$/)
		{
		my $query_file_gz = $query_file;
		$query_file =~ s/.gz//;
		system("$uncompress $query_file_gz");
		}
		
	push(@stats_vector, $query_file, "");

	# Count reads in query file
	my ($fa_total, $fa_unique) = (0, 0);
	if($fasta_format eq "r")	{	$fa_total=fastaTools::count_reads_fasta($query_file, "r"); $fa_unique="unknown";	}
	elsif($fasta_format eq "nr")	{	($fa_total, $fa_unique)=fastaTools::count_nr_fasta_reads($query_file, "nr");	}
	elsif($fasta_format eq "mnr")	{	($fa_total, $fa_unique)=fastaTools::count_mnr_fasta_reads($query_file, "mnr");	}

	push(@stats_vector, $fa_total, $fa_unique);

	# Loop over targets
	for(my $e=0; $e<=$#targets; $e++)
		{
		my @mismatches = @{$targets[$e]};
		my $target_path = shift(@mismatches);
		my @path_parts = split("/", $target_path);
		my $target_file=$path_parts[$#path_parts];

		# Loop over mismatches in target
		for(my $f=0; $f<=$#mismatches; $f++)
			{
			my $mm = $mismatches[$f];			
			my $outname = "${target_file}_"."${query_file}_"."${mm}mm".".pat";
			$outname =~ s/.fa//g;
		
			# Align query to target with current number of mismatches
			my ($countT, $countU, $complex, $cop_per_read, $countM)=misc::patman_align($query_path, $outname, $target_path, 3000000, $mm, $fasta_format);
			my $perc_total = ($countT/$fa_total)*100;
			my $perc_unique = ($countU/$fa_unique)*100;
			push(@stats_vector, "", $countT, $perc_total, $countU, $perc_unique, $complex, $cop_per_read, $countM);
			if(($zip eq "all") or ($zip eq "out"))	{	system("$compress $outname");	}
			}
		} # Loop over targets ends here
	my $outline = join(",", @stats_vector);
	print($stat "$outline\n");
	if(($zip eq "all") or ($zip eq "in"))	{	system("$compress $query_file");	}
	} # Loop over query files ends here

close($stat);
#close($log);
#close($wlog);
exit;

# end processing

########################################## Define local functions ##########################################

# end functions
