# DESCRIPTION: This script compares a list of transcripts assembled from RNAseq data ("data.gff") to an annotation file for the reference genome of the species of interest
# (also a *.gff file) and prints to an outfile all transcripts from the first file that don't fall within annotated features (= genes, mRNAs, rRNAs etc.)
# listed in the second file.
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

my $usage = "perl check_data_gff_against_annotated_gff.pl data.gff reference.gff outfile_name\n";
my $data = shift or die $usage;
my $ref = shift or die $usage;
my $outname = shift or die $usage;

if(!open(DATA, $data))       {       die "Couldn't open $data\n";        }
if(!open(REF, $ref))       {       die "Couldn't open $ref\n";        }

# Read the gff file into a matrix
my @raw_ref_matrix=();

# Format of reference.gff file:

#	0	seqid					NC_015762.1
#	1	method of annotation			Gnomon
#	2	annotation class			exon
#	3	start position				10000432
#	4	stop position				10000538
#	5	alignment score (not relevant)		.
#	6	strand					+
#	7	phase (of CDS:s)			.
#	8	attributes (various info)		ID=id8269;Parent=rna946;Dbxr... etc.

my $num_comment_lines=0;
while(<REF>)
	{
	my $line = text::trim($_);
	my @arr=split("\t", $line);
	foreach my $el (@arr)	{	$el = text::trim($el);	}
	if($line !~ /^#/)	{	push(@raw_ref_matrix, [($arr[0]), ($arr[3]), ($arr[4])]);	}
	else	{	$num_comment_lines++;	}
	}

my $num_ref_features=scalar(@raw_ref_matrix);
my $total_lines=$num_ref_features+$num_comment_lines;

# Format of @ref_matrix is now:

#	0	seqid		NC_015762.1
#	1	start		10000432
#	2	stop		10000538


# Sort @ref_matrix on seqid, start and stop
my @ref_matrix = sort { $a->[0] cmp $b->[0] || $a->[1] <=> $b->[1] || $a->[2] <=> $b->[2] } @raw_ref_matrix;

# Read the data.gff file into a matrix
my @raw_data_matrix=();

# Format of a line in data.gff:

#	0	seqid				NC_015762.1    
#	1	start				311
#	2	stop				2062
#	3	total abundance			551
#	4	positive strand abundance	548
#	5	negative strand abundance	3

while(<DATA>)
	{
	my $line=text::trim($_);
	my @arr=split("\t", $line);
	push(@raw_data_matrix, \@arr);
	}

# Sort the @raw_data_matrix on seqid, start and stop

my @data_matrix = sort {$a->[0] cmp $b->[0] || $a->[1] <=> $b->[1] || $a->[2] <=> $b->[2]} @raw_data_matrix;

# Open filehandle to an outfile that will hold the data transcripts that fall outside of annotated features
if(!open(OUT, ">>$outname"))	{	die "Couldn't create outfile una.gff\n";		}

# Compare the @data_matrix to the @ref_matrix

my $annotated=0;
my $novel=0;
my $total=scalar(@data_matrix);
my $index_ref=0;
my $stay_on_same_read="n";

# LOOP OVER DATA TRANSCRIPTS
DATA: for(my $k=0; $k<=$#data_matrix; $k++)
	{
	#if($stay_on_same_read eq "y")	{	print("(still) ");	}
	#print("Processing data transcript $k\t\t-\t");
	$stay_on_same_read="n";
	
	# IF DATA TRANSCRIPT ON EARLIER CHROMOSOME
	if($data_matrix[$k][0] lt $ref_matrix[$index_ref][0])
		{
		#print("Data transcript $k is on EARLIER chromosome than reference feature $index_ref\tnext data transcript!\n");

		# Add data transcript to the list of novel transcripts. Script goes to next data transcript. If last data transcript, loop ends.
		my @arr=@{$data_matrix[$k]};
		my $outline=join("\t", @arr);
		print(OUT "$outline\n");
		$novel++;
		}

	# IF DATA TRANSCRIPT ON SAME CHROMOSOME
	elsif($data_matrix[$k][0] eq $ref_matrix[$index_ref][0])
		{
		#print("Data transcript $k is on SAME chromosome as reference feature $index_ref\t");

		# IF DATA TRANSCRIPT BEFORE REF FEATURE
		if($data_matrix[$k][2] < $ref_matrix[$index_ref][1])
			{
			#print("and falls before reference feature $index_ref\tnext data transcript!\n");
			# Add data transcript to the list of novel transcripts. Script goes to next data transcript. If last data transcript, loop ends.
			my @arr=@{$data_matrix[$k]};
	                my $outline=join("\t", @arr);
        		print(OUT "$outline\n");
			$novel++;
			}
					
		# IF DATA TRANSCRIPT OVERLAPS REF FEATURE
		elsif(!($data_matrix[$k][1] > $ref_matrix[$index_ref][2]))
			{
			#print("and overlaps reference feature $index_ref\tnext data transcript!\n");
			$annotated++;				
			}

		# IF DATA TRANSCRIPT AFTER REFERENCE FEATURE
		elsif($data_matrix[$k][1] > $ref_matrix[$index_ref][2])
			{
			#print("and falls after reference feature $index_ref\t");
			# Unless we are on the last reference feature in the list, go to the next reference feature, but stay on the same data transcript
			if($index_ref != $#ref_matrix)
				{
				#print("next reference feature!\n");
				$index_ref++;
				$stay_on_same_read="y";
				$k--;				# Stay on the same data transcript
				# And now the script goes back to the top of the DATA loop,
				# where it now starts evaluating the same data transcript in relation to the new reference feature
				}
				
			# Else, if this is the last reference feature...
			elsif($index_ref == $#ref_matrix)
				{
				#print("which is the last reference feature\tnext data transcript!\n");
				# Add data transcript to list of novel transcripts. Script goes to next data transcripts. If last data transcripts, loop ends.
				my @arr=@{$data_matrix[$k]};
                             	my $outline=join("\t", @arr);
                           	print(OUT "$outline\n");
                        	$novel++;
				}		
			}
		}

	# DATA TRANSCRIPT ON LATER CHROMOSOME
	elsif($data_matrix[$k][0] gt $ref_matrix[$index_ref][0])
		{
		#print("Data transcript $k is on LATER chromosome than reference feature $index_gff\t");

		# Unless we are on the last reference feature in the list, go to the next reference feature, but stay on the same data transcript
		if($index_ref != $#ref_matrix)
			{
			#print("next reference feature!\n");
			$index_ref++;
                    	$stay_on_same_read="y";
			$k--;
                	# And now the script goes back to the top of the DATA loop,
                 	# where it now starts evaluating the same data transcript in relation to the new reference feature
			}
		# Else, if this is the last reference feature...
		elsif($index_ref == $#ref_matrix)
			{
			#print("which is the last reference feature\tnext data transcript!\n");
			# Add data transcript to the list of novel transcripts. Script goes to next data transcript. If last data transcript, loop ends.
                         my @arr=@{$data_matrix[$k]};
                         my $outline=join("\t", @arr);
                  	print(OUT "$outline\n");
           		$novel++;
			}
		}
	} # Here ends the loop over the @data_matrix (the transcripts assembled from the data)

print("\nFeatures in reference annotation:\t$num_ref_features\nComment lines in reference annotation:\t$num_comment_lines\nTotal lines:\t$total_lines\n\n");

if(!open(STAT, ">>stats_novel_transcripts.csv"))	{	die "Couldn't create statistics outfile\n";	}
print(STAT "Total_data_transcripts,Previously_annotated,Novel\n");
print(STAT "$total,$annotated,$novel\n");

close(OUT);
close(DATA);
close(REF);
close(STAT);
exit;

# end processing

########################################## Define local functions ##########################################

# end functions
