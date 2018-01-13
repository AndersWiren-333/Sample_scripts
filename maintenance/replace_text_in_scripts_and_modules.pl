# DESCRIPTION: This script searches all scripts and modules for a specific text string and optionally replaces it with another string.
# It is useful e.g. when renaming subroutines in modules (and you have to replace its name in all scripts and modules that use it).
#
# Example usage:	"perl replace_text_in_scripts_and_modules.pl $search_for $replace_with $test_or_sharp"
#
#	where	$search_for is the text to be searched for
#			$replace_with is the text to replace \$search_for with
#			$test_or_sharp  Set this to 'test' if you want to search a text but not replace it, or to 'sharp' if you want to do the replacement
#
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
sub process_dir;

# Define usage (error message to display if the user enters the wrong arguments at the command line)
my $usage = "Usage error for script ${0}. Correct usage: 'perl $0 \$search_for \$replace_with, \$test_or_sharp\n\nwhere".
"\t\$search_for is the text to be searched for\n".
"\t\$replace_with is the text to replace \$search_for with\n".
"\t\$test_or_sharp  Set this to 'test' if you want to search a text but not replace it, or to 'sharp' if you want to do the replacement\n\n";

# Accept input parameters
my @pars = @ARGV or die $usage;
foreach my $el (@pars)  {       $el = text::trim($el);  }
my $find = shift(@pars) or die $usage;
my $replace = shift(@pars) or die $usage;
my $test_or_sharp = shift(@pars) or die $usage;

# Go to scripts folder
chdir($scripts);
print "\nScripts\n----------------------------\n";
my @script_files=get_files_in_dir();
for(my $c=0; $c<=$#script_files; $c++)
	{
	my $file = $script_files[$c];
	#print "$file\n";
	process_file($file);
	}

# Go to modules folder
chdir($modules);
print "\nModules\n----------------------------\n";
my @module_files=get_files_in_dir();
for(my $c=0; $c<=$#module_files; $c++)
	{
	my $module = $module_files[$c];
	process_file($module);
	}

my $commit_message = "\n----------------------------------------\n".
"Suggested commit message:\n\nReplaced/renamed subroutine '$find' with/to '$replace' in all scripts and modules that use it.\n".
"Also replaced the name in 'Modules_and_functions.xlsx'\n\n\tNB! NB! Replacing in 'Modules_and_functions.xlsx' has to be done manually! NB!\n";	
if($test_or_sharp eq "sharp")	{	print "$commit_message\n";	}

#close($log);
#close($wlog);
exit;

# end processing

########################################## Define local functions ##########################################

sub get_files_in_dir
	{
	# List things in the current directory
	opendir(my $dh, ".");
	my @things=readdir($dh);
	my @files=();

	# Put files in a separate array
	for(my $cc=0; $cc<=$#things; $cc++)
		{
		my $thing = $things[$cc];
		if((-f $thing) && ($thing !~ /^\./))	{	push(@files, $thing);	}
		}
	closedir($dh);
	return(@files);
	}

sub process_file
	{
	# Loop over files
	my @pars = @_;
	my $infile = shift @pars;
	
	open(my $in, "<", $infile) or die "Script $0 couldn't open infile $infile\n";
	my $outname = "new_$infile";
	my @infile_cont = <$in>;
	close($in);
	my $infile_contents = join("", @infile_cont);
	
	print "Searched file $infile\t";
	if($infile_contents =~ /$find/)
		{
		print "found '$find'\t";
		if($test_or_sharp eq "sharp")
			{
			$infile_contents =~ s/$find/$replace/g;
			print "replaced it with '$replace'";
			open(my $out, ">", $outname) or die "Script $0 couldn't create outfile $outname\n";	
			print($out "$infile_contents");
			close($out);
			unlink($infile);
			move($outname, $infile);
			
			}
		}
	else	{	print "\t";	}	
	print "\n";
	
	#close($out);
	}

# end functions
