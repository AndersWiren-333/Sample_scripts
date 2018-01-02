# DESCRIPTION:
# This script takes 
#
#	a) a patman outfile (a file containing information about chromosome, start- and stop position, strand and abundances of sequence
#	   reads aligning to a reference genome).
#
#	b) an "una.gff" file containing the chromosome, start- and stop position, strand etc. of transcripts assembled from a specific
#	   RNAseq experiment, that do not overlap with any previously annotated features in the reference genome for the species of interest
#	   (= unannotated transcripts)
#
# The script assigns each read in the patman file to one of these unannotated transcripts, calculates the abundance of each transcript
# (= total number of reads mapping to that transcript) and writes the transcript ID and abundance to an outfile (one outfile for each
# patman input file).
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
my $usage = "perl novel_assign.pl una.gff genome_aligning_reads_sampleX.pat outfile_name.csv\n";
my $gff=shift or die $usage;
my $infile=shift or die $usage;
my $outname = shift or die $usage;
if(!open(GFF, $gff))	{	die "Couldn't open $gff\n";	}
if(!open(IN, $infile))  {       die "Couldn't open $infile\n";    }

my @raw_gff_matrix=();
my @raw_sample_matrix=();
my @out_matrix=();

# Read una.gff file into a matrix

# Format of a line in una.gff:

#	0	seqid (chromosome)		NC_015762.1
#	1	start				300
#	2	stop				600
#	3	total abundance			551
#	4	positive strand abundance	548
#	5	negative strand abundance	3

READ_GFF: while(<GFF>)
	{
	my $line=text::trim($_);
	my @arr = split("\t", $line);
	my $id=join("_", ($arr[0], $arr[1], $arr[2]));
	push(@raw_gff_matrix, [$arr[0],$arr[1],$arr[2],0,0,0,$id]);	# The zeros serve to reset the abundances of the gff features - those were the abundances of the feature when counting all samples,
	}								# but we want to count abundances for this sample only.
close(GFF);

# Sort @raw_gff_matrix
my @gff_matrix = sort {$a->[0] cmp $b->[0] || $a->[1] <=> $b->[1] || $a->[2] <=> $b->[2]} @raw_gff_matrix;

# Format of a line in @gff_matrix is now:

#       0       seqid (chromosome)              NC_015762.1
#       1       start                           300
#       2       stop                            600
#       3       total abundance                 551
#       4       positive strand abundance       548
#       5       negative strand abundance       3
#	6	transcript_id			NC_015762.1_300_600


# Read sample patman file into matrix

# Format of a line in the patman file:

#	0	seqid					NC_015762.1 Bombus terrestris linkage group LG B01, Bter_1.0, whole genome shotgun sequence
#	1	sequence-abundance			GTTCGCAGTACGCTCGCGATACATTCGTCGATCTCGAATCGATCGCGGAT-1
#	2	start					407
#	3	stop					456
#	4	strand					+
#	5	no. of mismatches in alignment		0

READ_SAMP: while(<IN>)
	{
	my $line=text::trim($_);
	my @arr = split("\t", $line);
	foreach my $el (@arr)	{	$el = text::trim($el);	}
	my @id_parts=split(" ", $arr[0]);
	my $id=$id_parts[0];
	my @seq_parts=split("-", $arr[1]);
	my $abund=$seq_parts[1];
	push(@raw_sample_matrix, [$id, $arr[2], $arr[3], $abund, $arr[4]]);	
	}
close(IN);

# Sort @raw_sample_matrix
my @sample_matrix = sort { $a->[0] cmp $b->[0] || $a->[1] <=> $b->[1] || $a->[2] <=> $b->[2] } @raw_sample_matrix;

# Format of a line in @sample_matrix is now:

#	0	seqid		NC_015762.1
#	1	start		407
#	2	stop		456
#	3	abundance	1
#	4	strand		+

# Assign reads to transcripts

my $index_ref=0;
my $stay_on_same_read="n";

# LOOP OVER PATMAN READS
DATA: for(my $k=0; $k<=$#sample_matrix; $k++)
        {
        #if($stay_on_same_read eq "y")  {       print("(still) ");      }
        #print("Processing read $k\t\t-\t");
        $stay_on_same_read="n";

        # IF READ ON EARLIER CHROMOSOME (than current transcript in gff matrix)
	print("$sample_matrix[$k][0]\t$gff_matrix[$index_ref][0]\n");
        if($sample_matrix[$k][0] lt $gff_matrix[$index_ref][0])
                {
                #print("read $k is on EARLIER chromosome than reference feature $index_ref\tnext read!\n");
                # Skip read (Script goes to next read). If last read, loop ends.
                }

        # IF READ ON SAME CHROMOSOME
        elsif($sample_matrix[$k][0] eq $gff_matrix[$index_ref][0])
                {
                #print("Read $k is on SAME chromosome as reference feature $index_ref\t");

                # IF READ BEFORE REF FEATURE
                if($sample_matrix[$k][2] < $gff_matrix[$index_ref][1])
                        {
                        #print("and falls before reference feature $index_ref\tnext read!\n");
                        # Skip read (Script goes to next read). If last data read, loop ends.
                        }

                # IF READ OVERLAPS REF FEATURE
                elsif(!($sample_matrix[$k][1] > $gff_matrix[$index_ref][2]))
                        {
                        #print("and overlaps reference feature $index_ref\tnext read!\n");
                        # Add the abundance of the current read to the abundance of the current gff feature
			$gff_matrix[$index_ref][3] = $gff_matrix[$index_ref][3]+$sample_matrix[$k][3];
			if($sample_matrix[$k][4] eq "+")	{	$gff_matrix[$index_ref][4] = $gff_matrix[$index_ref][4]+$sample_matrix[$k][3];	}
			if($sample_matrix[$k][4] eq "-")        {       $gff_matrix[$index_ref][5] = $gff_matrix[$index_ref][5]+$sample_matrix[$k][3];  }
                        }

                # IF READ AFTER REFERENCE FEATURE
                elsif($sample_matrix[$k][1] > $gff_matrix[$index_ref][2])
                        {
                        #print("and falls after reference feature $index_ref\t");
                        # Unless we are on the last reference feature in the list, go to the next reference feature, but stay on the same read
                        if($index_ref != $#gff_matrix)
                                {
                                #print("next gff feature!\n");
                                $index_ref++;
                                $stay_on_same_read="y";
                                $k--;                           # Stay on the same data transcript
                                # And now the script goes back to the top of the DATA loop,
                                # where it now starts evaluating the same read in relation to the new gff feature
                                }

                        # Else, if this is the last gff feature...
                        elsif($index_ref == $#gff_matrix)
                                {
                                #print("which is the last reference feature\tNow we are finished!\n");
                                # We have assigned the relevant reads to all gff features, so now we are finished! Exit the DATA loop.
				$k=$#sample_matrix;
                                }
                        }
                }

        # READ ON LATER CHROMOSOME
        elsif($sample_matrix[$k][0] gt $gff_matrix[$index_ref][0])
                {
                #print("Read $k is on LATER chromosome than gff feature $index_ref\t");

                # Unless we are on the last gff feature in the list, go to the next gff feature, but stay on the same read
                if($index_ref != $#gff_matrix)
                        {
                        #print("next gff feature!\n");
                        $index_ref++;
                        $stay_on_same_read="y";
                        $k--;
                        # And now the script goes back to the top of the DATA loop,
                        # where it now starts evaluating the same data read in relation to the new gff feature
                        }
                # Else, if this is the last gff feature...
                elsif($index_ref == $#gff_matrix)
                        {
                        #print("which is the last gff feature\tNow we are finished!\n");
             		# We have assigned the relevant reads to all gff features, so now we are finished! Exit the DATA loop.
               		$k=$#sample_matrix;
                        }
                }
        } # Here ends the loop over the @sample_matrix


# Print the data to an outfile
if(!open(OUT, ">>$outname"))	{	die "Couldn't create outfile";	}

GFF: for(my $m=0; $m<=$#gff_matrix; $m++)
	{
	my $seqid=$gff_matrix[$m][6];
	my $abundance=$gff_matrix[$m][3];
	print(OUT "u$seqid,$abundance\n");
	}

close(OUT);
exit;

# end processing

########################################## Define local functions ##########################################

# end functions
