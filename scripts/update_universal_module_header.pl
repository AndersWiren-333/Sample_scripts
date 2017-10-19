# DESCRIPTION: This script updates the universal perl script header of specified scripts using
# the perl template script (template_script.pl) and the update list (script_update_list.txt)
# end description

# Load libraries that this script depends on
use warnings;
use strict;
use File::Copy;

# Declare local functions
sub endings;

# Define usage (error message to display if the user enters the wrong arguments at the command line)
my $usage = "Syntax error for script ${0}. Correct usage: 'perl $0'\n";

chdir("../modules");
my $infiles = "module_update_list.txt";
open(my $template, "<", "template_module.pm") or die "Script $0 couldn't open template_module.pm\n";

# Read the template file in order to get the universal perl script header
my @template_cont = <$template>;
my $template_contents = join("", @template_cont);
close($template);
my $template_header="";
$template_contents=endings($template_contents);

if($template_contents =~ m/(#+ Universal perl module header.+end header)/s)    {       $template_header = $1;	}

# Read names of the modules that should be updated, one by one
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

# Loop over modules to be updated
for(my $c=0; $c<=$#update_list; $c++)
	{
	my $module = $update_list[$c];
	if($module =~ /\#/)	{	next;	}	# Exclude modules that should be excluded
	open(my $in, "<", $module) or die "Script $0 couldn't open infile $module\n";
	open(my $out, ">", "new_${module}") or die "Script $0 couldn't create outfile new_${module}\n";

	# Read the contents of the current module file
	my @module_cont = <$in>;
	my $module_contents = join("", @module_cont);
	$module_contents=endings($module_contents);

	# Divide the contents into parts
	my $module_top="";
	my $module_header="";
	my $module_functions="";

	if($module_contents =~ m/(package.+end\ssub\slist)/s)	{	$module_top=$1;	}
	if($module_contents =~ m/(#+ Universal perl module header.+#\send\sheader)/s)	{	$module_header=$1;	}
	if($module_contents =~ m/(#+ Functions.+#\send\sfunctions)/s)	{	$module_functions=$1;	}
	
	print($out "$module_top\n\n");
	print($out "$template_header\n\n");
	print($out "$module_functions\n");

	close($in);
	close($out);

	unlink($module);
	move("new_${module}", $module);

	} # Loop over modules ends here

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
