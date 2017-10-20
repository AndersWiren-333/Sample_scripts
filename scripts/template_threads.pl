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

use threads;

# Declare local functions (if any)
sub do_something_locally;
sub call_script;

# Define usage (error message to display if the user enters the wrong arguments at the command line)
my $usage = "Usage error for script ${0}. Correct usage: 'perl $0 \$infile \$infile, \$maxlines'\n\nwhere".
"\t\$infile is the file that should be used as input\n".
"\t\$maxlines is the maximum number of lines the infile may have without being split into parts\n\n";
my @pars = @ARGV or die $usage;
foreach my $el (@pars)  {       $el = text::trim($el);  }
my $infile = shift(@pars) or die $usage;
my $maxlines = shift(@pars) or die $usage;

my ($partnames_ref, $number_of_parts)=text::split_file($infile, $maxlines);
my @partnames = @{$partnames_ref};

# Have $number_of_parts threads filter the partfiles in parallel
my $from=1;
my $to=5;
my $num_threads=3;
my @threads=(1..$num_threads);	# Create thread id vector
print("\nStarting $num_threads processes... \n");

#foreach(@threads)	{	$_ = threads->create(\&do_something_locally);	}	# I guess I could call another script here... ?
#foreach(@threads)      {       $_ = threads->create('do_something_locally', $from, $to);  }
foreach(@threads)      {       $_ = threads->create("call_script", $from, $to, "utfil_");  }
foreach(@threads)	{	$_->join();	}

print("Finished!\n");


#close($in);
#close($out);
#close($stat);
#close($log);
#close($wlog);
exit;

# end processing

########################################## Define local functions ##########################################

sub do_something_locally
	{
	my $id = threads->tid();	# Create a new thread id (?)
	my $from=shift;
	my $to=shift;
	open(my $out, ">", "outfile_${id}.txt") or die "Couldn't create thread outfile\n";
	for(my $i=$from; $i<=$to; $i++)	{	print($out "$i\n");	}
	close($out);
	}

sub call_script
	{
	my $id = threads->tid();
        my $from=shift;
        my $to=shift;
	my $outnamebase=shift;
	my $outname = "$outnamebase"."${id}.txt";
	system("perl $scripts/zzz.pl $from $to $outname");
	}

# end functions
