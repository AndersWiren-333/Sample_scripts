# DESCRIPTION: This is a reference collection of perl commands for working with databases
# (specifically SQLite)
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
#use diagnostics;	# Turn this on if you need more verbose warnings and error messages
use FindBin;
use DBI;	# This is needed for working with databases

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

# Declare local functions (if any)

# Define usage (error message to display if the user enters the wrong arguments at the command line)
my $usage = "Usage error for script ${0}. Correct usage: 'perl $0 \$infile \$namefile, \$outname'\n\nwhere".
"\t\$infile is the file that should be used as input\n".
"\t\$namefile is a text file listing names of the files, one on each line, that should be used as input files\n".
"\t\$outname is the desired name of the outfile that will be produced\n\n";
# my @pars = @ARGV or die $usage;
# foreach my $el (@pars)  {       $el = text::trim($el);  }
# my $infile = shift(@pars) or die $usage;
# my $namefile = shift(@pars) or die $usage;
# my @rest_of_pars = @pars or die $usage;

# Set variables
my $database = "databas.db";
my $dsn="dbi:SQLite:dbname=$database";	# Set the data source name
my $user="";
my $pass="";

# Open database (in SQLite this is automatically creates a database if it doesn't already exist)
my $dbhandle = DBI->connect($dsn, $user, $pass);

# Create table
$dbhandle->do("create table Alver(first_name text, last_name text, born_age text, born_year int)");

my $fname='Bilbo';
my $lname='Baggins';
my $age='Third Age';
my $year=1343;


# Insert data into table
$dbhandle->do('insert into Alver (first_name, last_name, born_age, born_year) values(?,?,?,?)', undef, $fname, $lname, $age, $year);

my $newyear=312;
my $name='Galadriel';

# Update fields
$dbhandle->do('update Alver set born_year = ? where first_name = ?', undef, $newyear, $name);

# Select something
my $sql='select first_name, last_name from Alver where born_age = ?';	# Prepare statement
my $sthandle = $dbhandle->prepare($sql);		# Open a statement handle
$sthandle->execute('First Age');			# Execute the select query_string

# Capture the output of the query to an array for each row
while(my @row = $sthandle->fetchrow_array)
	{
	print "$row[0] $row[1]\n";
	}

print "\n\n";
	
# Capture the output of the query to a hash
$sthandle->execute('Third Age');
while(my $row = $sthandle->fetchrow_hashref)
	{
	print "$row->{first_name} $row->{last_name}\n";
	}
	
$dbhandle->disconnect;


#open(my $in, "<", $infile) or die "Script $0 couldn't open infile $infile\n";
#open(my $out, ">", $outname) or die "Script $0 couldn't create outfile $outname\n";
#open(my $stat, ">", "stats.csv") or die "Script $0 couldn't create statistics outfile\n";





# ALL OTHER CODE HERE

#close($in);
#close($out);
#close($stat);
#close($log);
#close($wlog);
exit;

# end processing

########################################## Define local functions ##########################################

# end functions
