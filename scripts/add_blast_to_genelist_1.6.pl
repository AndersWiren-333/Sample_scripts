# DESCRIPTION: This script adds blast information to a list of differentially expressed genes (in csv format)
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
my $usage = "Usage error for script ${0}. Correct usage: 'perl $0 \$genelist.csv \$blast_hit_file.txt \$list_of_blasted_sequences \$outfile_name'\n\nwhere".
"\t\$genelist.csv is the genelist that should be added to\n".
"\t\$blast_hit_file.txt is a textfile listing all blasthits (one on each line) for all blasted sequences\n".
"\t\$list_of_blasted_sequences is a textfile listing (one on each line) the Feature_IDs of all the sequences that have been blasted\n".
"\t\$outfile_name.csv is the desired name of the outfile\n\n";

# Declare variables and filehandles
my $genelist=shift or die $usage;
my $blastfile=shift or die $usage;
my $blasted_file = shift or die $usage;
my $outfile=shift or die $usage;


if(!open(LIST, $genelist))	{	die "Couldn't open $genelist\n";	}
if(!open(BLAST, $blastfile))      {       die "Couldn't open $blastfile\n";        }
if(!open(ORIG, $blasted_file))	{	die "Couldn't open list of blasted sequences $blasted_file\n";	}
if(!open(OUT, ">>$outfile"))        {       die "Couldn't open outfile\n";  }

my @genelist_headers=();
my @genelist=();
my $blast_start;
my @ofc_indices=();

# Read list of blasted sequences into vector
my @orig_blasted=();

while(<ORIG>)
	{
	my $line = trim($_);
	push(@orig_blasted, $line);
	}
close(ORIG);


# Read genelist into matrix

my $gene_counter=0;
GENES: while(<LIST>)
	{
	my $line = trim($_);
	my @arr = split(",", $line);
	for(my $i=0; $i<=$#arr; $i++)	{	$arr[$i]=trim($arr[$i]);	}
	
# Format of arr is now:

#	0	Feature_ID	NC_015770.1_12249827_12250301
#	1-x	Sample data	225.1
#	x+1	length		475
#	x+2	OFC_comp_1	1.23
#	x+3	OFC_comp_2	0.89
#	x+4	seqid		NC_015770.1
#	x+5	source		RNAseq
#	x+6	type		novel_transcript
#	x+7	start		12249827
#	x+8	end		12250301
#	x+9	score	
#	x+10	strand		+/-
#	x+11	phase	
#	x+12	BLAST_hit_01	no hit
#	x+13	BLAST_hit_02	no hit
#	x+14	BLAST_hit_03	no hit
#	x+15	ID_01		rna21647
#	x+16	Name_01		XM_012320077.1
#	x+17	Product_01	cuticle collagen 2C-like

# etc...

	if($gene_counter==0)	# If this is the header line...
		{
		@genelist_headers=@arr;		# Set headers
		($blast_start) = grep { $genelist_headers[$_] eq "BLAST_hit_01" } 0..$#genelist_headers;
		(@ofc_indices) = grep { $genelist_headers[$_] =~ /OFC/ } 0..$#genelist_headers;			
		$gene_counter++;
		next GENES;
		}
	else		# If it is any other row...
		{
		push(@genelist, \@arr);
		$gene_counter++;
		}
	}
close(LIST);

@genelist = sort { $a->[0] cmp $b->[0] } @genelist;


# Read Blast-file into matrix (but only top three hits for each gene)

my @blast_hits=();
BLAST: while(<BLAST>)
	{
	my $line = trim($_);
	my @arr=split("\t", $line);
	for(my $i=0; $i<=$#arr; $i++)    {       $arr[$i]=trim($arr[$i]);        }
	$arr[1] =~ s/,//g;
	$arr[0] =~ s/u//;
	push(@blast_hits, \@arr);
	}


@blast_hits = sort { $a->[0] cmp $b->[0] || $b->[3] <=> $a->[3] || $b->[6] <=> $a->[6] || $b->[5] <=> $a->[5] || $a->[2] <=> $b->[2] } @blast_hits;

# The format of a line in @blast_hits is now:

#	0	Query seqid		(first sorting criterion)
#	1	Subject title		
#	2	E-value			(fifth sorting criterion)
#	3	Bitscore		(second sorting criterion)
#	4	Raw score
#	5	Alignment length	(fourth sorting criterion)
#	6	% identity		(third sorting criterion)

# Remove all entries from the blast matrix that contain the word "uncharacterized" (those entries are uninformative
# and only drown out the information we really need)


# Loop over sorted list of blast hits.
my @new_blast_hits=();
REDUCE: for(my $i=0; $i<=$#blast_hits; $i++)
	{
	if($blast_hits[$i][1] =~ /uncharacterized/)	{	next REDUCE;	}
	else	{	push(@new_blast_hits, $blast_hits[$i]);	}
	}
@blast_hits=@new_blast_hits;
@new_blast_hits=();

# Reformat the blast hit matrix so that each gene has one entry, and the format of that entry is this:

#       0       Gene_id
#       1       Name of blast-hit 1
#       2       Name of blast-hit 2
#       3       Name of blast-hit 3

my @new_blast_matrix=();	# Make new matrix to hold the formatted list
my $comp_id="";
my @outarr=();	# Initiate a small array to hold a gene id and that genes' top three blast hits


BLAST_MATRIX: for(my $i=0; $i<=$#blast_hits; $i++)
	{
	my @arr=@{$blast_hits[$i]};
	my $id=$arr[0];	
	
	if($i==0)		# If this is the first line of the hit matrix...
		{
		$comp_id=$id;
		@outarr=($arr[0], $arr[1]);		# Start outline with Feature_ID and hit name (since this is the first (=top) hit for that gene)	
		}
	else			# If this is any other line...
		{
		if($id eq $comp_id)		# If we are still on the same gene...
			{
			push(@outarr, $arr[1]);		# Add the name of the current hit to the growing hit list for the current gene
			if($i==$#blast_hits)			# If this is the last line of the hit matrix...
				{
				my @newarr=@outarr[0..3];		# ...add the top three blast hits for the current gene to the new, reformatted blast matrix
				push(@new_blast_matrix, \@newarr);
				}
			}
		else				# If we are on a new gene...
			{
			my @newarr=@outarr[0..3];		# add the top three blast hits for the previous gene to the new, reformatted blast matrix
			push(@new_blast_matrix, \@newarr);
			$comp_id=$id;			# ...start growing the hit list for the new gene.
                	@outarr=($arr[0], $arr[1]);
			if($i==$#blast_hits)
				{
				push(@new_blast_matrix, \@newarr);
				}
			}

		}
		

	}

@new_blast_matrix = sort { $a->[0] cmp $b->[0] } @new_blast_matrix;

if(!open(MAT, ">>matrix.csv")) {	die "Couldn't create matrix.csv";	}	
for(my $ii=0; $ii<=$#new_blast_matrix; $ii++)
	{
	my @line = @{$new_blast_matrix[$ii]};
	my $outline = join(",", @line);
	print(MAT "$outline\n");
	}
close(MAT);

# Loop over genes in the genelist, and if the BLAST_hit_01 field says "do", add the top three blast hits to the
# record for that gene.

print("Entries in genelist: " . scalar(@genelist) . "\n");
print("Entries in blast_matrix: " . scalar(@new_blast_matrix) . "\n");

# Set new column names for blast hit columns
my $out_headers=join(",", @genelist_headers);
print(OUT "$out_headers\n");

my @out_matrix2=();

# Loop over genelist
GENELIST: for(my $m=0; $m<=$#genelist; $m++)
	{
	my @gene=@{$genelist[$m]};
	my $seqid=$gene[0];

	my $hit="n";
	my $blasted="n";

	# Loop over genes in reformatted blast matrix (check if current gene has a blast hit)
	for(my $n=0; $n<=$#new_blast_matrix; $n++)
		{
		my $blast_id=$new_blast_matrix[$n][0];

		# If the genelist gene id is the same as the current gene id in the blast matrix...			
		if($seqid eq $blast_id)
			{
			$hit="y";
			$gene[$blast_start] = $new_blast_matrix[$n][1];		# Add first blast hit to first blast hit column of genelist
			$gene[1+$blast_start] = $new_blast_matrix[$n][2];	# Add second blast hit to second blast hit column of genelist
			$gene[2+$blast_start] = $new_blast_matrix[$n][3];	# Add third blast hit to third blast hit column of genelist
			}
		} # Loop over blast hit matrix ends here


	# Loop over entries in @orig_blasted vector
	for(my $o=0; $o<=$#orig_blasted; $o++)
		{
		my $blasted_id = $orig_blasted[$o];
		if($seqid eq $blasted_id)	{	$blasted="y";	}
		}

	if(($hit eq "n") && ($blasted eq "y"))	{	$gene[$blast_start] = "no relevant hit";	}
	if(($hit eq "n") && ($blasted eq "n"))  {       $gene[$blast_start] = "not blasted";        }
	
	push(@out_matrix2, \@gene);
	}

my $ofc_1=$ofc_indices[0];
my @out_matrix3 = sort { $b->[$ofc_1] <=> $a->[$ofc_1] } @out_matrix2;

for(my $s=0; $s<=$#out_matrix3; $s++)
	{
	my @inline = @{$out_matrix3[$s]};
	my $outline = join(",", @inline);
	print(OUT "$outline\n");
	}

unlink("matrix.csv");
close(ERR);
close(OUT);
close(BLAST);
exit;

# end processing

########################################## Define local functions ##########################################

# end functions
