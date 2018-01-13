# DESCRIPTION: Use this script if you want to apply an update of a file in this repository to a corresponding file
# (with the same name) in one or more different repositories. The script uses the file "commit_messages.txt" for
# information about which file to update in which repositories, whether to add, update or delete the file and which
# commit message to use.
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
$scripts =~ s/\/maintenance/\/scripts/;
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
require "$modules/normalise.pm";
require "$modules/listTools.pm";

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

# Define usage (error message to display if the user enters the wrong arguments at the command line)
my $usage = "Usage error for script ${0}. Correct usage: 'perl $0'\n\n";

my $infile="";
my @repos=();

print "\n\n";

# Read commit message
my @commit_file=misc::read_list("../commit_messages.txt");
my $commit_str = join("\n", @commit_file);

# Get the commit message
my $message="";
my $mode="";
if($commit_str =~ /--message(.+)--end_message/s)	{	$message = $1;	}
$message =~ s/^\n//;
$message =~ s/\n$//;

# Get name of infile and repositories to be updated
for(my $c=0; $c<=$#commit_file; $c++)
	{
	my $line = $commit_file[$c];
	if($line =~ /file=/)
		{
		(my $part1, $infile) = split("=", $line);
		$infile = text::trim($infile);
		}
	elsif($line =~ /repo=/)
		{
		unless($line =~ /^#/)
			{
			my ($part2, $repo) = split("=", $line);
			$repo=text::trim($repo);
			push(@repos, $repo);
			}
		}
	elsif($line =~ /mode=/)
		{
		(my $part1, $mode) = split("=", $line);
		$mode = text::trim($mode);
		}
	}
	
# If the infile is a specific (perl, for the time being) subroutine
my $module="";
my $sub="";
if($infile =~ /::/)
	{
	($module, $sub) = split("::", $infile);
	my $modpath = "${modules}/${module}.pm";
	
	# Read correct subroutine from the infile
	print "Opening $modpath to read subroutine $sub\n";
	open(my $modfile, "<", $modpath) or die "Script $0 couldn't open infile $module\n";
	my @modlines = <$modfile>;
	my $modstring = join("", @modlines);
	close($modfile);
	
	# Get the relevant subroutine
	my $substring="";
	if($modstring =~ /(.+)(sub ${sub}\s.+ # end $sub\s)(.+)/s)	{	$substring = $2; print "Got substring from infile module\n";	}
	else	{	print "Subroutine isn't present in module. \$mode is $mode\n";	}
	
	# Loop over replacement repos
	for(my $d=0; $d<=$#repos; $d++)
		{
		# Read correct file
		my $repository = $repos[$d];
		print "Processing repository $repository\n";
		chdir("${repository}/modules"); # Changing to repo module directory
		my $modpath2 = "${repository}/modules/${module}.pm";
		
		print "Opening $modpath2 to read subroutine $sub\n";
		open(my $modfile2, "<", $modpath2) or die "Script $0 couldn't open infile $modpath2\n";
		my @modlines2 = <$modfile2>;
		my $modstring2 = join("", @modlines2);
		close($modfile2);
		
		# Replace subroutine
		my $first_part2="";
		my $substring2="";
		my $last_part2="";
		my $newmod="${module}.pm";
		print "Deleting old file $modpath2\n";
		unlink($newmod);
		
		open(my $outmod, ">>", $newmod) or die "Script $0 couldn't create outfile $newmod\n";
		
		# If the subroutine should be deleted
		if($mode eq "delete")
			{
			if($modstring2 =~ /(.+)(sub ${sub}\s.+ # end $sub\s)(.+)/s)
				{
				$first_part2 = $1;
				$substring2 = $2;
				$last_part2 = $3;
				}
			print($outmod "$first_part2");
			print($outmod "$last_part2");
			}
		
		# If this is an update to an already existing subroutine
		elsif($modstring2 =~ /(.+)(sub ${sub}\s.+ # end $sub\s)(.+)/s)
			{
			$first_part2 = $1;
			$substring2 = $2;
			$last_part2 = $3;
			
			print($outmod "$first_part2");
			print($outmod "$substring\n\n");
			print($outmod "$last_part2");
			}
			
		# If it is a new subroutine (= it is at the end of its module)
		elsif($modstring2 !~ /(sub ${sub}\s.+ # end $sub\s)/s)
			{
			if($modstring2 =~ /(.+)(return\(1\).+# end functions)/s)
				{
				$first_part2 = $1;
				$last_part2 = $2;
				}
			
			print($outmod "$first_part2");
			print($outmod "$substring\n\n");
			print($outmod "$last_part2");
			}

		close($outmod);

		chdir(".."); # Changing to main repo directory
		system("git add modules/${newmod}");
		system("git commit -m \"${message}\"");
		}
	}
	
# If the infile is a whole file
else
	{
	my $source_repo = $modules;
	$source_repo =~ s/\/modules//;
	my $source_file = "${source_repo}/${infile}";
	
	# Loop over replacement repos
	for(my $e=0; $e<=$#repos; $e++)
		{
		my $dest_repo = $repos[$e];
		chdir($dest_repo); # Changing to main repo directory
		my $dest_file = "${dest_repo}/${infile}";
		
		# If a new file should be added
		if($mode eq "add")
			{
			copy($source_file, $dest_file);
			system("git add $infile");
			system("git commit -m \"${message}\"");
			}
		
		# If the file should be deleted
		elsif($mode eq "delete")
			{
			unlink($dest_file);
			system("git rm $infile");
			system("git commit -m \"${message}\"");
			}
		
		# If the file should be updated
		elsif($mode eq "update")
			{
			unlink($dest_file);
			copy($source_file, $dest_file);
			system("git add $infile");
			system("git commit -m \"${message}\"");
			}		
		}
	}

#close($log);
#close($wlog);
exit;

# end processing

########################################## Define local functions ##########################################

# end functions
