# DESCRIPTION: This script takes a list of differentially expressed genes and a number of patman output files (genome alignment information for the sequence reads
# that the differential gene expression ananlysis is based on). For each gene in the list, the script calculates - for each nucleotide position along
# the gene - the sequencing depth at that position. For eaxample, if three sequencing reads overlap position 342 in the gene, the depth for that
# position is three. A separate depth is calculated for the positive and negative strand of the gene. The script writes the information to a csv outfile,
# which should be used as input for the script make_presence_plots.r

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

# Declare variables and filehandles
my $usage = "perl get_presence_plot_data.pl genelist.csv namefile num_genes";	# namefile is a list of patman outfiles, one on each line. num_genes is the number of (top DE) genes to get data for.
my $genelist = shift or die $usage;
my $namefile = shift or die $usage;
my $num_genes_to_print = shift or die $usage;
my $limit=1+$num_genes_to_print;
if(!open(GENE, $genelist))	{	die "Couldn't open $genelist\n";	}
if(!open(NAM, $namefile))      {       die "Couldn't open $namefile\n";        }
my $outname = "presence_data_".$genelist;
if(!open(OUT, ">>$outname"))      {       die "Couldn't open $outname\n";        }
my @out_matrix=();


# Read genelist into matrix (the top ten (or $limit) genes)
my @genelist_matrix=();
my $genecounter=0;
my $seqid_field=0;
my $start_field=0;
my $stop_field=0;

RGENES: while(<GENE>)
	{
	$genecounter++;
	if($genecounter<=$limit)
		{
		my $line = text::trim($_);
		my @arr=split(",", $line);
		if($genecounter==1)		# If this is the first line, find indices for seqid, start and stop
			{
			($seqid_field) = grep { $arr[$_] eq "seqid" } 0..$#arr;
			($start_field) = grep { $arr[$_] eq "start" } 0..$#arr;
			($stop_field) = grep { $arr[$_] eq "end" } 0..$#arr;
			}
		# If not, put the genes chromosome, start and stop position in a new matrix that becomes a line of @genelist_matrix
		else	{
			$arr[$seqid_field] =~ s/u//;
			my @newarr=($arr[$seqid_field], $arr[$start_field], $arr[$stop_field]);
			push(@genelist_matrix, \@newarr);
			}
		}
	else	{	last RGENES;	}
	}

# Format of @genelist_matrix is now:

#	0	seqid (chromosome)
#	1	start
#	2	stop

# Read the namefile into an array of patman filenames

my @filenames=();
while(<NAM>)
	{
	my $file=text::trim($_);
	push(@filenames, $file);
	}

# Start looping over samples (patman output files)
foreach my $el (@filenames)
	{
	my $filename=$el;
	my $sampname=$filename;
	$sampname =~ s/genome_//;
	$sampname =~ s/_1mm.pat//;
	if(!open(PAT, $filename))	{	die "Couldn't open $filename\n";	}

	# Read the patman file into a matrix
	my @pat_matrix=();
	while(<PAT>)
		{
		my $line=text::trim($_);
		my @arr=split("\t", $line);
		my @ids=split(" ", $arr[0]);		# HAR KAN DET MOJLIGEN BLI PROBLEM! (eftersom forsta elementet i raden (seqid) i just de har patmanfilerna inte innehaller nagra blanksteg) 
		$arr[0]=$ids[0];
		my @seqparts=split("-", $arr[1]);
		$arr[1]=$seqparts[1];
		push(@pat_matrix, \@arr);
		}
	
	# Format of @pat_matrix is now:

	#	0	seqid
	#	1	abundance
	#	2	start
	#	3	stop
	#	4	strand
	#	5	mismatches



	# For each gene in the genelist.... 

	for(my $i=0; $i<=$#genelist_matrix; $i++)
		{
		my $generank=1+$i;
		my $geneid=$genelist_matrix[$i][0];
		my $genestart=$genelist_matrix[$i][1];
		my $genestop=$genelist_matrix[$i][2];
		my $genename=$geneid."_".$genestart."_".$genestop;
		my $length=$genestop-$genestart+1;
		my @geneheader=($genestart..$genestop);
		my @geneheader2=(1..$length);
		my @gene_pos = (0) x $length;
		my @gene_neg = (0) x $length;

		# Check each read in @pat_matrix....

		READS: for(my $j=0; $j<=$#pat_matrix; $j++)
			{
			my $readid=$pat_matrix[$j][0];
			my $readstart=$pat_matrix[$j][2];
			my $readstop=$pat_matrix[$j][3];
			my $readstrand=$pat_matrix[$j][4];
			my $readabund=$pat_matrix[$j][1];

			if($readid eq $geneid)
				{
				if($readstart >= $genestart)
					{
					if($readstop <= $genestop)
						{
						my ($startind) = grep { $geneheader[$_] == $readstart } 0..$#geneheader;
						my ($stopind) = grep { $geneheader[$_] == $readstop } 0..$#geneheader;
						if($readstrand eq "+")	
							{
							for(my $k=$startind; $k<=$stopind; $k++)	{	$gene_pos[$k]=$gene_pos[$k]+$readabund;	}
							}
						elsif($readstrand eq "-")
							{
                                                        for(my $k=$startind; $k<=$stopind; $k++)        {       $gene_neg[$k]=$gene_neg[$k]-$readabund; }
                                                        }
						}
					}
				}
			}
		my $outheader=join(",", @geneheader2);
		my $outline_pos=join(",", @gene_pos);
		my $outline_neg=join(",", @gene_neg);
		push(@out_matrix, [$generank, $genename, $sampname, "pos", $outline_pos]);
		push(@out_matrix, [$generank, $genename, $sampname, "neg", $outline_neg]);
		}
	}

# Format of out_matrix is now:

#	0	gene rank
#	1	gene name
#	2	sample name	
#	3	strand
#	4	data

# Sort @out_matrix on gene rank, sample name and strand

@out_matrix = sort { $a->[0] <=> $b->[0] || $a->[2] cmp $b->[2] || $b->[3] cmp $a->[3] } @out_matrix;

# Print @out_matrix to a csv outfile (for later import to R)

for(my $m=0; $m<=$#out_matrix; $m++)
	{
	my @arr=@{$out_matrix[$m]};
	my $outline=join(",", @arr);
	print(OUT "$outline\n");
	}

close(GENE);
close(NAM);
close(OUT);
exit;

# end processing

########################################## Define local functions ##########################################

# end functions
