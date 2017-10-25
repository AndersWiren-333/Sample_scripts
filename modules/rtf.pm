package rtf;

# openRtf($filename)
# closeRtf($filename)
# stringRtf($string, $colour, $bold, $italic, $underline, $font, $fontsize)
# printRtf($string)
# end sub list

########################################################## Universal perl module header ##########################################################

# Load libraries that this module depends on
use warnings;
use strict;
use Sys::Hostname;
use File::Copy;
use Bio::SeqIO;
use Cwd;
use threads;
use diagnostics;
use FindBin;

# Set paths to scripts and modules. Setting explicit paths to the scripts and modules in this specific repository (rather than adding paths to @INC, PERLLIB and PATH on
# your local system) avoids the risk of scripts calling the wrong scripts/modules if you have other repositories on your system that happen to have some script- and module names
# in common with this repository.
my $thisfile = (__FILE__);
my $modules = "";
if($thisfile =~ m/^(.+)\//)	{	$modules = $1;	}
my $scripts = $modules;
$scripts =~ s/modules/scripts/;

# If this script/module is intended to be used outside the folder structure of the parent repository (e.g. a wrapper script to be started from
# another part of your system), set the absolute path to repository scripts and modules (that this cript may depend on) here (and comment out
# the lines for seeting paths above). Otherwise, keep the lines below commented out.
#my $modules = "path/to/modules/folder";
#my $scripts = "path/to/scripts/folder";

# Load home-made modules (aka "libraries", aka "packages") that this module depends on
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

# end header

########################################################## Functions ##########################################################

# Starts an rtf file. After this you can use stringRtf to write text to the file, and last you have to use closeRtf to close the file.
sub openRtf
	{
	my $usage = "Syntax error for sub openRtf. Correct usage:\n\t'rtf::openRtf(\$filename)'\n";
	my $filename = $_[0] or die $usage;
	$filename = "${filename}.rtf";
	if(!open(RTF, ">>$filename"))	{	die "Couldn't create rtf file $filename with sub openRtf\n";	}
	
	# Open file and set font table
	print(RTF "{\\rtf1\\ansi\\deff0{\\fonttbl\n");	# Start rtf prolog
	print(RTF "{\\f0 Courier New;}\n");	# Set font 0 of the font table
	print(RTF "{\\f1 Calibri;}\n");		# Set font 1 of the font table
	print(RTF "{\\f2 Times New Roman;}\n");		# Set font 2 of the font table
	print(RTF "}\n\n");		# End font table
	
	# Set colour table
	print(RTF "{\\colortbl;");	# Start colour table and set colour 0 (black?)
	print(RTF "\\red0\\green0\\blue0;");	# Set colour 1 (black)
	print(RTF "\\red153\\green0\\blue204;");	# Set colour 2 (purple)
	print(RTF "\\red0\\green0\\blue255;");	# Set colour 3 (dark blue)
	print(RTF "\\red0\\green176\\blue240;");	# Set colour 4 (blue)
	print(RTF "\\red0\\green255\\blue255;");	# Set colour 5 (turquoise)
	print(RTF "\\red0\\green153\\blue0;");	# Set colour 6 (dark green)
	print(RTF "\\red0\\green255\\blue0;");	# Set colour 7 (light green)
	print(RTF "\\red255\\green255\\blue0;");	# Set colour 8 (yellow)
	print(RTF "\\red255\\green192\\blue0;");	# Set colour 9 (orange)
	print(RTF "\\red255\\green0\\blue0;");	# Set colour 10 (red)
	print(RTF "\\red165\\green0\\blue33;");	# Set colour 11 (dark red)
	
	print(RTF "}\n\n");		# End colour table
	
	# Print a document formatting header
	print(RTF "{\\header\\pard\\qr\\plain\\f0\\par}\n\n");
	
	print(RTF "{\\pard ");	# Start a paragraph
	}

# Ends an rtf file. You have to use this as the last thing you do/print to an rtf file, or it will not be readable.
sub closeRtf
	{
	my $usage = "Syntax error for sub closeRtf. Correct usage:\n\t'rtf::closeRtf(\$filename)'\n";
	print(RTF "\\par}\n\n}");
	}
	
# Creates a string formatted in rtf format, using the given text string, colour (black, purple, dark blue, blue, tuquoise, dark green, light green, yellow, orange, red, dark red), given
# boldness (yes or no), italicision (yes/no), underlining (yes/no), font (times, courier, calibri) and font size. You have to first open the file with openRtf and when you have created a
# formatted string with this function, print it to the open rtf file using printRtf. When you are finished, don't forget to close the file with closeRtf)
sub stringRtf
	{
	my $usage = "Syntax error for sub stringRtf. Correct usage:\n\t'rtf::stringRtf(\$string, \$colour, \$bold, \$italic, \$underline, \$font, \$fontsize)'\n";
	my $string = $_[0] or die $usage;
	my $col = $_[1] or die $usage;
	my $bold = $_[2] or die $usage;
	my $italic = $_[3] or die $usage;
	my $underline = $_[4] or die $usage;
	my $font = $_[5] or die $usage;
	my $fontsize = $_[6] or die $usage;
	
	# Convert fontsize to rtf "twips" units (= 2* normal Microsoft fontsize)
	$fontsize = 2*$fontsize;
	$fontsize = "\\fs$fontsize ";
	
	# Replace normal newlines and tabs with rtf newlines and tabs
	$string =~ s/\n/\\line /g;
	$string =~ s/\t/\\tab /g;
	
	# Set colour
	if($col eq "black")	{	$col = "\\cf1 "; }
	elsif($col eq "purple")	{	$col = "\\cf2 "; }
	elsif($col eq "dark blue")	{	$col = "\\cf3 "; }
	elsif($col eq "blue")	{	$col = "\\cf4 "; }
	elsif($col eq "turquoise")	{	$col = "\\cf5 "; }
	elsif($col eq "dark green")	{	$col = "\\cf6 "; }
	elsif($col eq "light green")	{	$col = "\\cf7 "; }
	elsif($col eq "yellow")	{	$col = "\\cf8 "; }
	elsif($col eq "orange")	{	$col = "\\cf9 "; }
	elsif($col eq "red")	{	$col = "\\cf10 "; }
	elsif($col eq "dark red")	{	$col = "\\cf11 "; }
	
	# Set boldness, italicision and underlining
	if($bold eq "no")	{	$bold = "";	}
	elsif($bold eq "yes")	{	$bold = "\\b ";	}
	if($italic eq "no")	{	$italic = "";	}
	elsif($italic eq "yes")	{	$bold = "\\i ";	}
	if($underline eq "no")	{	$underline = "";	}
	elsif($underline eq "yes")	{	$underline = "\\ul ";	}
	
	# Set font
	if($font eq "courier")	{	$font = "\\f0 ";	}
	elsif($font eq "calibri")	{	$font = "\\f1 ";	}
	elsif($font eq "times")	{	$font = "\\f2 ";	}
	
	my $new_string = ("$font"."$fontsize"."$col"."$bold"."$italic"."$underline"."$string");
	#print(RTF "$font"."$fontsize"."$col"."$bold"."$italic"."$underline"."$string");
	}

# Writes a pre-formatted rtf-string to an open rtf file. First open the file with openRtf, then create the pre-formatted rtf-string with stringRtf, then print
# it to the rtf file using this function and finally don't forget to close the rtf file with closeRtf when you are finished.
sub printRft
	{
	my $usage = "Syntax error for sub printRtf. Correct usage:\n\t'rtf::printRtf(\$string)'\n";
	my $string = $_[0] or die $usage;
	print(RTF "$string");
	}

return(1);

# end functions
