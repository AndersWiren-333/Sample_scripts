package dnaTools;

# sub get_revcomp($sequence)
# sub get_comp($sequence)
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

# Constructs the reverse complement of a DNA or RNA sequence, with ambiguity codes.
sub get_revcomp
        {
	my $usage = "Syntax error for sub get_revcomp. Correct usage: 'dnaTools::get_revcomp(\$sequence);'\n";
        my $seq = $_[0] or die $usage;
        my $last=(length($seq))-1;
        my @seqarr = split("", $seq);
        my @new_arr=();

	# Loop backwards over nucleotides in input sequence
        for(my $cc=$last; $cc>=0; $cc--)
                {
                my $nuc=$seqarr[$cc];
                my $new_nuc="";
                if($nuc eq "A") {	$new_nuc = "T"; }
                elsif($nuc eq "C") {	$new_nuc = "G"; }
                elsif($nuc eq "G") {	$new_nuc = "C"; }
                elsif($nuc eq "T") {	$new_nuc = "A"; }
		elsif($nuc eq "U") {	$new_nuc = "A"; }
                elsif($nuc eq "N") {	$new_nuc = "N"; }
                elsif($nuc eq "-") {	$new_nuc = "-"; }

                elsif($nuc eq "M")	{	$new_nuc = "K"; }
                elsif($nuc eq "R")	{	$new_nuc = "Y"; }
                elsif($nuc eq "W")	{	$new_nuc = "W"; }
                elsif($nuc eq "S")	{	$new_nuc = "S"; }
                elsif($nuc eq "Y")	{	$new_nuc = "R"; }
                elsif($nuc eq "K")	{	$new_nuc = "M"; }

		elsif($nuc eq "B")	{	$new_nuc = "V";	}
		elsif($nuc eq "D")	{	$new_nuc = "H";	}
		elsif($nuc eq "H")	{	$new_nuc = "D";	}
		elsif($nuc eq "V")	{	$new_nuc = "B";	}
				
                push(@new_arr, $new_nuc);
                }
        my $new_seq = join("", @new_arr);
        return($new_seq);
        }


# Constructs the complement of a DNA or RNA sequence, with ambiguity codes.
# Parameters: $sequence
sub get_comp
        {
        my $usage = "Syntax error for sub get_comp. Correct usage: 'dnaTools::get_comp(\$sequence);'\n";
        my $seq = $_[0] or die $usage;
        my @seqarr = split("", $seq);
        my @new_arr=();

  	# Loop over nucleotides in input sequence
        for(my $cc=0; $cc<=$#seqarr; $cc++)
                {
                my $nuc=$seqarr[$cc];
                my $new_nuc="";
                if($nuc eq "A") {       $new_nuc = "T"; }
                elsif($nuc eq "C") {    $new_nuc = "G"; }
                elsif($nuc eq "G") {    $new_nuc = "C"; }
                elsif($nuc eq "T") {    $new_nuc = "A"; }
             	elsif($nuc eq "U") {    $new_nuc = "A"; }
                elsif($nuc eq "N") {    $new_nuc = "N"; }
                elsif($nuc eq "-") {    $new_nuc = "-"; }

                elsif($nuc eq "M")      {       $new_nuc = "K"; }
                elsif($nuc eq "R")      {       $new_nuc = "Y"; }
                elsif($nuc eq "W")      {       $new_nuc = "W"; }
                elsif($nuc eq "S")      {       $new_nuc = "S"; }
                elsif($nuc eq "Y")      {       $new_nuc = "R"; }
                elsif($nuc eq "K")      {       $new_nuc = "M"; }

             	elsif($nuc eq "B")      {       $new_nuc = "V"; }
           	elsif($nuc eq "D")      {       $new_nuc = "H"; }
            	elsif($nuc eq "H")      {       $new_nuc = "D"; }
              	elsif($nuc eq "V")      {       $new_nuc = "B"; }

                push(@new_arr, $new_nuc);
                }
        my $new_seq = join("", @new_arr);
        return($new_seq);
        }

return(1);

# end functions
