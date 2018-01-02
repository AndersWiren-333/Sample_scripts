# DESCRIPTION:
# This script 
#	a) finds information about genes in an expression matrix (where gene names are in the form NCxxxxxx.x_start-pos_stop-pos) from the reference genome gff file
#	b) finds the sequences of unannotated and uncharacterized genes/transcripts from the reference genome fasta file, and prints them to a new output fasta file
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
my $usage = "Usage: perl get_info_from_gff.pl gene_expression_matrix.csv annotation_file.gff sample_name_1 sample_name_2 ...";
my $infile = shift or die $usage;
my $annotation = shift or die $usage;
my @samples=@ARGV;
if(!open(IN, $infile))	{	die "Couldn't open $infile";	}
if(!open(OUT, ">>an1_$infile"))	{	die "Couldn't create outfile";	}	# File to write matches to, along with information
if(!open(AN, $annotation))	{	die "Couldn't open $annotation";	}


# Read the annotation (gff) file into a matrix and format it in a better way

my @annotation=();

while(<AN>)
	{
	my $line = text::trim($_);
	if($line =~ /^#/)	{	next;	}	# Remove the comment rows at the beginning of the gff file
	
	my @arr = split("\t", $line);
	for(my $a=0; $a<=$#arr; $a++)	{	$arr[$a]=text::trim($arr[$a]);	}	

	my $new_id = join("_", @arr[0,3,4]);    # Make a new ID for the gene (Chromosome_start_stop)
	unshift(@arr, $new_id);                # Add the new ID to the beginning of the array

	# The format of @arr is now:
	# newID         seqid   source  type    start   end     score   strand  phase   attributes
	# 0             1       2       3       4       5       6       7       8       9


	# Split up the attributes and keep only the relevant parts
	my @attributes=();
	my @attribute_parts = split(";", $arr[9]);

	my @things=("","","");

	foreach my $attribute_part (@attribute_parts)
		{
		if($attribute_part =~ /(ID=.+)/)        {       my @mini=split("=", $1); $things[0]=$mini[1];   }
		if($attribute_part =~ /(Name=.+)/)      {       my @mini=split("=", $1); $things[1]=$mini[1];      }
		if($attribute_part =~ /(product=.+)/)      {       my @mini=split("=", $1); $things[2]=$mini[1];      }
		}
	my $zero_things=0;
	foreach my $thing (@things)     {       if(($thing eq "") or (!defined $thing) )        {       $zero_things++; }       }
	if($zero_things != 3)   {       push(@attributes, @things);       }

	@arr = (@arr[0..8], , "", "", "", @attributes);

	# The format of @arr is now:
	# newID         seqid   source  type    start   end     score   strand  phase   BLAST_hit_01	BLAST_hit_02	BLAST_hit_03	ID_01	Name_01	 Product_01    ID_02   Name_02	Product_02	etc.
	# 0             1       2       3       4       5       6       7       8       9		10		11	      	12      13	 14		15	16	17

	push(@annotation, [@arr]);		# Add the new line to the genome matrix
	}
close(AN);

# Print headers to outfile
my @headers= ("Feature_ID",@samples,
		"length","seqid","source","type","start","end","score","strand","phase","BLAST_hit_01","BLAST_hit_02","BLAST_hit_03",
		"ID_01","Name_01","Product_01","ID_02","Name_02","Product_02","ID_03","Name_03","Product_03","ID_04","Name_04","Product_04","ID_05","Name_05","Product_05",
		"ID_06","Name_06","Product_06","ID_07","Name_07","Product_07","ID_08","Name_08","Product_08","ID_09","Name_09","Product_09","ID_10","Name_10","Product_10",);
my ($source_field) = grep { $headers[$_] eq "source"  } 0..$#headers;
my ($blast_start) = grep { $headers[$_] eq "BLAST_hit_01" } 0..$#headers;
my ($ID_02) = grep { $headers[$_] eq "ID_02" } 0..$#headers;
my ($seqid_field) = grep { $headers[$_] eq "seqid" } 0..$#headers;
my $out_headers=join(",", @headers);
print(OUT "$out_headers\n");

# Read the gene expression matrix line by line and match the gene identifier to those in the annotation file (gff file)
# If there is a match, print all the info to an outfile
while(<IN>)
	{
	my $inline=text::trim($_);
	my @arr=split(",", $inline);	
	for(my $b=0; $b<=$#arr; $b++)	{	$arr[$b] = text::trim($arr[$b]);	}

	my $inid="";
	my $instart="";
	my $instop="";

	if($arr[0] =~ /(.+)(_)(\d+)(_)(\d+)/)
		{	
		$inid = $1;
		$instart = $3;
		$instop = $5;
		}

	my $length=($instop-$instart)+1;
	my $newinid=join("_", $inid, $instart, $instop);
	
	# Match gene in matrix to reference genome annotation...
	my $hit_counter=0;

	# For every line in the gff file...
	my @hits=();					# @hits will hold zero or more hits for this specific gene (it is a local list) 
	for(my $i=0; $i<=$#annotation; $i++)
		{
		if($newinid eq $annotation[$i][0])
			{
			push(@hits, $annotation[$i]);	# $genome[$i] is an array reference
			$hit_counter++;
			}
		}
	
	my @outarr=();		# @outarr will hold all the data to be printed to the outfile for this specific gene

	# If there are several hits for each gene, put the information from each hit on the same line, but without
	# unnecessary repetition, just the ID, Name and Product fields.
	if(scalar(@hits) == 1)
		{
		my @newinfo=@{$hits[0]};
		shift(@newinfo);
		@outarr = (@arr, $length, @newinfo);
		}
	elsif(scalar(@hits == 0))
		{
		@outarr = (@arr, $length, $inid, "RNAseq", "novel_transcript", $instart, $instop);
		}
	else
		{
		my @newinfo=@{$hits[0]};
                shift(@newinfo);
		@outarr = (@arr, $length, @newinfo);
		for(my $j=1; $j<scalar(@hits); $j++)	
			{
			my @hitarr=@{$hits[$j]};
			push(@outarr, @hitarr[9..$#hitarr]);
			}
		}

	# Format of @outarr:
	#	0	Feature_ID
	# 	1-18	Raw data
	#	19	length
	#	20	seqid
	#	21	source
	#	22	type
	#	23	start
	#	24	end
	#	25	score
	#	26	strand
	#	27	phase
	#	28	BLAST_hit_01
	#	29	BLAST_hit_02
	#	30	BLAST_hit_03
	#	31	ID_01
	#	32	Name_01
	#	33	Product_01
	#	34	ID_02
	#	35	Name_02
	#	36	Product_02
	#	etc...

	# Check if the current gene is unannotated and needs to be blasted
	my $novel=0;
	my $undefined=1;
	my $unchar=0;
	if($outarr[$source_field] eq "RNAseq")	{	$novel=1;	}
	else
		{
		UNCHAR: for(my $i=$ID_02; $i<=$#outarr; $i++)
			{
			if($outarr[$i] ne "")	{	$undefined=0;	}
			if($outarr[$i] =~ /uncharacterized/)	{	 $unchar=1;	}
			}
		}
	if(($novel==1) or ($undefined==1) or ($unchar==1))	# If the gene needs blasting...
		{
		$outarr[$blast_start]="do";				# Label it (in the column that will later hold the first BLST hit)
		}

	# Print the resulting line to the outfile, in csv format
	my $outline=join(",", @outarr);
	print(OUT "$outline\n");
	}

close(IN);
close(OUT);
exit;

# end processing

########################################## Define local functions ##########################################

# end functions
