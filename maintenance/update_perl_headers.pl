# DESCRIPTION: This script updates the universal perl header of scripts and modules conatining the flag "perl_script_update"
# or "perl_module_update" in their header. Note that they must contain the universal script or module header to begin with.
# end description

#!/usr/bin/perl

# Load libraries that this script depends on
use warnings;
use strict;
use File::Copy;

# Declare local functions
sub trim;
sub get_files_in_dir;
sub process_scripts;
sub process_modules;

# Define usage (error message to display if the user enters the wrong arguments at the command line)
my $usage = "Syntax error for script ${0}. Correct usage: 'perl $0\n";
open(my $script_template, "<", "template_script.pl") or die "Script $0 couldn't open script_template_script.pl\n";
my $script_flag = "perl_script_update";
my $module_flag = "perl_module_update";

# Read the template script to get the universal perl script header
my @script_template_cont = <$script_template>;
my $script_template_contents = join("", @script_template_cont);
close($script_template);
my $script_template_header="";
if($script_template_contents =~ m/(#+ Universal perl script header.+end header)/s)    {       $script_template_header = $1;	}

# Read the template module to get the universal perl script header
open(my $module_template, "<", "template_module.pm") or die "Script $0 couldn't open template_module.pm\n";
my @module_template_cont = <$module_template>;
my $module_template_contents = join("", @module_template_cont);
close($module_template);
my $module_template_header="";
if($module_template_contents =~ m/(#+ Universal perl module header.+end header)/s)    {       $module_template_header = $1;	}

print "\n";

# Go to scripts folder and process scripts
chdir("../scripts");
my @scripts = get_files_in_dir();
process_scripts(\@scripts);
print "\n";

# Go to maintainance folder and process maintainance scripts
chdir("../maintainance");
my @maintainance_scripts = get_files_in_dir();
process_scripts(\@maintainance_scripts);
print "\n";

# Go to module folder and process modules
chdir("../modules");
my @modules = get_files_in_dir();
process_modules(\@modules);
print "\n";

#close($log);
#close($wlog);
exit;
# end processing

########################################## Define local functions ##########################################

sub trim
	{
	# Set error messages and accept input parameters
	my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
	my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$textstring)'\n\nwhere".
	"\t\$textstring is the string of text that should be trimmed\n\n";
	my $string = shift or 0;

	$string =~ s/\r+//;
	$string =~ s/\n+//;
	$string =~ s/^\s+//m;
	$string =~ s/\s+$//m;
	return($string);
	}

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
	
sub process_scripts
	{
	# Set error messages
	my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
	my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$scriptnames_array_ref)'\n\nwhere".
	"\t\$scriptnames_array_ref is a reference to the array containing the names of scripts to be processed\n\n";
	
	# Accept input parameters
	my @pars = @_ or die $usage;
	foreach my $el (@pars)  {       $el = trim($el);  }
	my $update_list_ref = shift @pars or die $usage;	
	my @update_list = @{$update_list_ref};
	
	# Loop over scripts to be updated
	for(my $c=0; $c<=$#update_list; $c++)
		{
		my $script = $update_list[$c];
		open(my $in, "<", $script) or die "Script $0 couldn't open infile $script\n";

		# Read the contents of the current script file
		my @script_cont = <$in>;
		my $script_contents = join("", @script_cont);
		close($in);
		unless($script_contents =~ /$script_flag/)	{	next;	}
		if(($script eq "update_perl_headers.pl") or ($script eq "template_script.pl"))	{	next;	}
		
		# Divide the contents into parts
		my $script_description="";
		my $script_header="";
		my $script_processing="";
		my $script_functions="";

		if($script_contents =~ m/(# DESCRIPTION.+end description)/s)	{	$script_description=$1;	}
		if($script_contents =~ m/(#+ Universal perl script header.+end header)/s)	{	$script_header=$1;	}
		if($script_contents =~ m/(#+ Processing.+end processing)/s)	{	$script_processing=$1;	}
		if($script_contents =~ m/(#+ Define local functions.+end functions)/s)	{	$script_functions=$1;	}

		open(my $out, ">", "new_${script}") or die "Script $0 couldn't create outfile new_${script}\n";
		print($out "$script_description\n\n");
		print($out "$script_template_header\n\n");
		print($out "$script_processing\n\n");
		print($out "$script_functions\n");
		close($out);

		unlink($script);
		move("new_${script}", $script);
		print "Processed script $script\n";
		} # Loop over scripts ends here	
	}
	
sub process_modules
	{
	# Set error messages
	my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
	my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$modulenames_array_ref)'\n\nwhere".
	"\t\$modulenames_array_ref is a reference to the array containing the names of modules to be processed\n\n";
	
	# Accept input parameters
	my @pars = @_ or die $usage;
	foreach my $el (@pars)  {       $el = trim($el);  }
	my $update_list_ref = shift @pars or die $usage;	
	my @update_list = @{$update_list_ref};
	
	# Loop over modules to be updated
	for(my $c=0; $c<=$#update_list; $c++)
		{
		my $module = $update_list[$c];
		open(my $in, "<", $module) or die "Script $0 couldn't open infile $module\n";

		# Read the contents of the current module file
		my @module_cont = <$in>;
		close($in);
		my $module_contents = join("", @module_cont);
		#print "$module\n";
		unless($module_contents =~ /$module_flag/)	{	next; print "\tskipped $module\n";	}
		if($module eq "template_module.pm")	{	next; print "\tskipped $module\n";	}

		# Divide the contents into parts
		my $module_top="";
		my $module_header="";
		my $module_functions="";

		if($module_contents =~ m/(package.+end\ssub\slist)/s)	{	$module_top=$1;	}
		if($module_contents =~ m/(#+ Universal perl module header.+#\send\sheader)/s)	{	$module_header=$1;	}
		if($module_contents =~ m/(#+ Functions.+#\send\sfunctions)/s)	{	$module_functions=$1;	}
		
		open(my $out, ">", "new_${module}") or die "Script $0 couldn't create outfile new_${module}\n";
		print($out "$module_top\n\n");
		print($out "$module_template_header\n\n");
		print($out "$module_functions\n");
		close($out);

		unlink($module);
		move("new_${module}", $module);
		print "Processed module $module\n";
		} # Loop over modules ends here	
	}
	
# end functions
