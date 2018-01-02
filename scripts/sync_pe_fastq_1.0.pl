# DESCRIPTION: This script synchronizes paired end fastq files (it discards unpaired reads)
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

# Declare local functions (if any)
sub read_fastq;
sub testprint_matrix;
sub print_matrix;

# Declare variables and filehandles
my $usage = "perl sync_fastq.pl file_1.fastq file2.fastq";
my $file_01 = shift or die $usage;
my $file_02 = shift or die $usage;
my $newname_01 = "new_$file_01";
my $newname_02 = "new_$file_02";
my $sing_01 = "singletons_$file_01";
my $sing_02 = "singletons_$file_02";

# Read files into matrices
my @one = read_fastq($file_01);
my @two = read_fastq($file_02);

@one = sort { $a->[0] cmp $b->[0] } @one;
@two = sort { $a->[0] cmp $b->[0] } @two;

my @singletons_01=();
my @singletons_02=();


# Compare matrices side by side, step by step

my @long_matrix=();
my @short_matrix=();

if(scalar(@one) >= scalar(@two))
	{
	@long_matrix = @one;
	@short_matrix = @two;
	}

elsif(scalar(@two) > scalar(@one))
	{
	@long_matrix = @two;
	@short_matrix = @one;
	}



# Compare matrices

my $counter_01=0;
my $counter_02=0;
my @new_01=();
my @new_02=();


MATLINES: while($counter_01<=$#one)
	{
	my $id_01="";
	my $id_02="";
	
	if($one[$counter_01][0] =~ /(@\w+:\d{2}:\w+:\d:\d+:\d+:\d+) (\d:\w:\d:\w{6})/)	{	$id_01=$1;	}
	if($two[$counter_02][0] =~ /(@\w+:\d{2}:\w+:\d:\d+:\d+:\d+) (\d:\w:\d:\w{6})/)	{	$id_02=$1;	}

	if($id_01 ne $id_02)
		{
		# ... if id 1 is earlier in alphabetic order than id 2...
		while($id_01 lt $id_02)
			{
			push(@singletons_01, $one[$counter_01]);		# Add that sequence to the singletons list for matrix @one
			$counter_01++;
			if($one[$counter_01][0] =~ /(@\w+:\d{2}:\w+:\d:\d+:\d+:\d+) (\d:\w:\d:\w{6})/)   {       $id_01=$1;      }	# Update $id_01 to the next sequence id
			}

		# ... else, if id 2 is earlier in alphabetic order than id 1...
		while($id_01 gt $id_02)
			{	
			push(@singletons_02, $two[$counter_02]);		# Add that sequence to the singletons list for matrix @two
			$counter_02++;
			if($two[$counter_02][0] =~ /(@\w+:\d{2}:\w+:\d:\d+:\d+:\d+) (\d:\w:\d:\w{6})/)   {       $id_02=$1;      }       # Update $id_02 to the next sequence id
			}
		}
	elsif($id_01 eq $id_02)
		{
		push(@new_01, $one[$counter_01]);
		push(@new_02, $two[$counter_02]);
		$counter_01++;
		$counter_02++;
		}


	}

my $ett = scalar(@new_01);
my $tva = scalar(@new_02);

print("\n$ett\t$tva\n\n");


print_matrix(\@new_01, $newname_01);
print_matrix(\@new_02, $newname_02);
print_matrix(\@singletons_01, $sing_01);
print_matrix(\@singletons_02, $sing_02);

exit;

# end processing

########################################## Define local functions ##########################################

sub read_fastq
	{
	my $file = shift;
	if(!open(FIL, $file))	{	die "Couldn't open $file_01\n";	}
	my @matrix=();
	my @temp=();
	my $counter=0;
	while(<FIL>)
		{
		$counter++;
		my $line = text::trim($_);
		push(@temp, $line);
		if($counter==4)				# If this is the last line of a sequence record...
			{
			push(@matrix, [@temp]);		# Add the four lines of the record as an element (a line of) to @one matrix
			@temp=();					# Free the @temp array
			$counter=0;				# Reset the counter
			}
		}
	return(@matrix);
	}
	
sub testprint_matrix
	{
	my $matrix_ref = shift;
	my @matrix = @{$matrix_ref};
	for(my $c=0; $c<=$#matrix; $c++)
		{
		my @arr=@{$matrix[$c]};
		foreach my $el (@arr)	{	print("$el\n");	}
		#print("\n---------------------\n\n");
		}
	}
	
sub print_matrix
	{
	my $matrix_ref = shift;
	my $outname = shift;
	my @matrix = @{$matrix_ref};
	if(!open(OUT, ">>$outname"))	{	die "Couldn't create outfile $outname\n";	}
	for(my $c=0; $c<=$#matrix; $c++)
		{
		my @arr=@{$matrix[$c]};
		foreach my $el (@arr)	{	print(OUT "$el\n");	}
		#print("\n---------------------\n\n");
		}
	close(OUT);
	}

# end functions
