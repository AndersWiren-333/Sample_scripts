# DESCRIPTION: This script updates the universal perl script header of specified scripts using
# the perl template script (template_script.pl) and the update list (script_update_list.txt)
# end description

#!/usr/bin/perl

# Load libraries that this script depends on
use warnings;
use strict;
use File::Copy;

# Declare local functions
sub endings;

# Define usage (error message to display if the user enters the wrong arguments at the command line)
my $usage = "Syntax error for script ${0}. Correct usage: 'perl $0\n";

my $infiles = "script_update_list.txt";
open(my $template, "<", "template_script.pl") or die "Script $0 couldn't open template_script.pl\n";

# Read the template file in order to get the universal perl script header
my @template_cont = <$template>;
my $template_contents = join("", @template_cont);
close($template);
my $template_header="";

$template_contents=endings($template_contents);
if($template_contents =~ m/(#+ Universal perl script header.+end header)/s)    {       $template_header = $1;	}

# Read names of the scripts that should be updated, one by one
my @update_list=();

open(my $filelist, "<", $infiles) or die "Script $0 couldn't open $infiles\n";

while(my $file = <$filelist>)
		{
		$file =~ s/\r+//;
		$file =~ s/\n+//;
		$file =~ s/\s+//m;
		$file =~ s/\s+$//m;
		push(@update_list, $file);
		}	

# Loop over scripts to be updated
for(my $c=0; $c<=$#update_list; $c++)
	{
	my $script = $update_list[$c];
	if($script =~ /\#/)	{	next;	}	# Exclude scripts that should be excluded
	open(my $in, "<", $script) or die "Script $0 couldn't open infile $script\n";
	open(my $out, ">", "new_${script}") or die "Script $0 couldn't create outfile new_${script}\n";

	# Read the contents of the current script file
	my @script_cont = <$in>;
	my $script_contents = join("", @script_cont);
	$script_contents=endings($script_contents);

	# Divide the contents into parts
	my $script_description="";
	my $script_header="";
	my $script_processing="";
	my $script_functions="";

	if($script_contents =~ m/(# DESCRIPTION.+end description)/s)	{	$script_description=$1;	}
	if($script_contents =~ m/(#+ Universal perl script header.+end header)/s)	{	$script_header=$1;	}
	if($script_contents =~ m/(#+ Processing.+end processing)/s)	{	$script_processing=$1;	}
	if($script_contents =~ m/(#+ Define local functions.+end functions)/s)	{	$script_functions=$1;	}

	print($out "$script_description\n\n");
	print($out "$template_header\n\n");
	print($out "$script_processing\n\n");
	print($out "$script_functions\n");

	close($in);
	close($out);

	unlink($script);
	move("new_${script}", $script);

	} # Loop over scripts ends here

#close($log);
#close($wlog);
exit;
# end processing

########################################## Define local functions ##########################################

sub endings
	{
	my $string = shift;
	#$string =~ s/\r//g;
	#$string =~ s/\n/\r\n/g;
	return($string);
	}

# end functions
