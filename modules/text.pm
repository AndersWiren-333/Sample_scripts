package text;

# sub trim($string)
# sub linebreaks($file, $desired_linebreaks)
# sub parse_path($path)
# sub read_sample_list($sample_names)
# sub split_file($infile.txt, $lines_per_piece)
# sub concat_files($outname, $file1, $file2, $file3)
# sub remove_spaces_from_csv($infile.csv)
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

# Removes whitespace from the ends of a given string, as well as all sorts of linebreaks from anywhere in the string
# Parameters: $string
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

# Replaces all linebreaks in a file with the type of linebreaks ('win' for Windows, 'lin' for Linux/Unix) you specify
# Parameters: $infile, $breakstyle (breakstyle = win or lin)
sub linebreaks
	{
        # Set error messages and accept input parameters
        my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
        my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$infile, \$desired_linebreaks)'\n\nwhere".
        "\t\$textstring is the string of text that should be trimmed\n".
	"\t\$desired_linebreaks can be either 'win' for windows or 'lin' for linux/unix\n\n";

	my $file = $_[0] or die $usage;
	my $osy = $_[1] or die $usage;
	$file = text::trim($file);
	$osy = text::trim($osy);
	open(my $in, "<", $file) or die "Subroutine $subname couldn't open infile $file\n";
	my $outname = ("$osy"."_"."$file");
	open(my $out, ">>", $outname) or die "Subroutine $subname couldn't create outfile $outname\n";

	my $end="";
	if($osy eq "win")	{	$end="\r\n";	}
	elsif($osy eq "lin")	{	$end="\n";	}
	else	{	die "Invalid breakstyle specified. Start again.\n";	}

	# Process file
	while(my $line = <$in>)
		{
		$line = text::trim($line);
		print($out "${line}${end}");
		}
	close($in);
	close($out);
	return($outname);
	}

# Gets the directory name, filename, file basename and file suffix of a path
sub parse_path
	{
        my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
        my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$path)'\n\nwhere".
        "\t\$path is a path\n\n";
        my @pars = @_ or die $usage;
        foreach my $el (@pars)  {       $el = text::trim($el);  }
        my $path = shift @pars or die $usage;

	my @parts = split("/", $path);
	my $file = pop(@parts);
	my $dir = join("/", @parts);
	my @fileparts = split(/\./, $file);
	my $suffix = pop(@fileparts);
	my $basename = join(".", @fileparts);
	
	my @outarr = ($dir, $file, $basename, $suffix);
        return(@outarr);
	}

# Reads a list of sample names into a matrix of sample group names. The sample name file must have the format 'samplename,group', e.g.
#
#	sample1,groupA
#	sample2,groupA
#	sample3,groupB
#	sample4,groupB
#
# The matrix that is returned has the format 'group,sample_index1,sample_index2,sample_index3', e.g.
#
#	groupA,0,1
#	groupB,2,3
sub read_sample_list
	{
	# Set error messages and accept input parameters
        my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
        my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$sample_list.txt)'\n\nwhere".
        "\t\$sample_list.txt is a plain text file with, on each line, the name of a sample followed by a comma and the name of the group to which that sample belongs, i.e.\n\n".
        "\t\tsample1,groupA\n\t\tsample2,groupA\n\t\tsample3,groupB\n\t\tsample4,groupB\n\n";
        my @pars = @_ or die $usage;
        foreach my $el (@pars)  {       $el = text::trim($el);  }
        my $namefile = shift @pars or die $usage;

	my $inmatrix_ref=(misc::file_to_matrix($namefile, "comma", "n"))[0];
	my @inmatrix=@{$inmatrix_ref};

	# Loop over samples in the @inmatrix while adding a numeric index to each sample and building a list of group names
	my $sample_index=0;
	for(my $c=0; $c<=$#inmatrix; $c++)
		{
		$inmatrix[$c][2]=$sample_index;
		$sample_index++;
		}

	# Sort the matrix on group name
	@inmatrix = sort { $a->[1] cmp $b->[1] || $a->[2] <=> $b->[2] } @inmatrix;

	# Build a new matrix based on group names
	my $comp_group = $inmatrix[0][1];
	my @temp_outline = ($inmatrix[0][1], $inmatrix[0][2]); 
	my @outmatrix=();
		
	for(my $d=1; $d<=$#inmatrix; $d++)
		{
		if($inmatrix[$d][1] eq $comp_group)
			{
			push(@temp_outline, $inmatrix[$d][2]);
			if($d==$#inmatrix)	{	push(@outmatrix, [@temp_outline]);	}
			}
		else
			{
			push(@outmatrix, [@temp_outline]);
			$comp_group = $inmatrix[$d][1];
			@temp_outline = ($inmatrix[$d][1], $inmatrix[$d][2]);
			if($d==$#inmatrix)	{	push(@outmatrix, [@temp_outline]);	}
			}
		}
	return(@outmatrix);
	}

# Splits a given data file into a number of pieces with a given number of lines in each, and returns a reference to an array holding the part names, and a scalar corresponding to the
# number of parts created. The last piece may have a smaller number. NB! If you want to split a fasta or fastq file you should not use this script.
# NB! If you want to split a fasta or fastq file you should not use this script.
sub split_file
	{
	# Set error messages and accept input parameters
        my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
        my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$infile.txt, \$lines_per_piece)'\n\nwhere".
        "\t\$infile.txt is the file that should be split\n".
        "\t\$lines_per_piece is the maximum number of lines in each outfile\n\n";
        my @pars = @_ or die $usage;
        foreach my $el (@pars)  {       $el = text::trim($el);  }
        my $infile = shift @pars or die $usage;      # NB! This will not work if the argument is the number 0 (because it will then be interpreted as false). In that case you need to use 'shift(@pars) or 0', but it may lead to proble$
	my $lines_per_file = shift @pars or die $usage;	

	my $partnum=1;
        open(my $in, "<", $infile) or die "Subroutine $subname (called by script '${calling_script}', line ${calling_line}) couldn't open infile $infile\n";

	my @partnames=();
	my $outname = "$infile"."_part_${partnum}";
	push(@partnames, $outname);
	open(our $out, ">", $outname) or die "Subroutine $subname (called by script '${calling_script}', line ${calling_line}) couldn't create outfile $partnum\n";
	my $linenum=1;

	while(my $line = <$in>)
		{
		$line = text::trim($line);
		if($linenum <= $lines_per_file)	{	print($out "$line\n");	}	
		else
			{
			close($out);
			$partnum++;
			$linenum=1;
			$outname = "$infile"."_part_${partnum}";
			push(@partnames, $outname);
			open(our $out, ">", $outname) or die "Subroutine $subname (called by script '${calling_script}', line ${calling_line}) couldn't create outfile $partnum\n";
			print($out "$line\n");
			}
		$linenum++;
		}
        close($in);
        close($out);
	my @results=(\@partnames, $partnum);
        return(@results);
	}

# Concatenates a series of data files.
sub concat_files
	{
        # Set error messages and accept input parameters
        my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
        my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$outname, \$file1, \$file2, \$file3 etc...)'\n\nwhere".
        "\t\$outname is the desired name of the resulting concatnated file\n".
        "\t\$file1, \$file2 etc are the files to be concatenated\n\n";
        my @pars = @_ or die $usage;
        foreach my $el (@pars)  {       $el = text::trim($el);  }
	my $outname = shift(@pars) or die $usage;

	open(my $out, ">>", $outname) or die "Subroutine $subname (called by script '${calling_script}', line ${calling_line}) couldn't create outfile $outname\n";

	for(my $cc=0; $cc<=$#pars; $cc++)
		{
		my $file = $pars[$cc];
		open(my $in, "<", $file) or die "Subroutine $subname (called by script '${calling_script}', line ${calling_line}) couldn't open infile $file\n";
		while(my $line = <$in>)	{	print($out "$line");	}
		close($in);
		unlink($file);	# This operation doesn't work for some reason (unlink is a built in function and usually works in this sort of context)
		}
	close($out);
	}

# Removes unwanted blank spaces from csv files
sub remove_spaces_from_csv
	{
        # Set error messages and accept input parameters
        my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
        my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$infile.csv)'\n\nwhere".
        "\t\$infile.csv is the infile to be cleaned\n\n";
        my @pars = @_ or die $usage;
        foreach my $el (@pars)  {       $el = text::trim($el);  }
        my $infile = shift @pars or die $usage;      # NB! This will not work if the argument is the number 0 (because it will then be interpreted as false). In that case you need to use 'shift(@pars) or 0', but it may lead to proble$
	my $outname = "t_".$infile;

        open(my $in, "<", $infile) or die "Subroutine $subname (called by script '${calling_script}', line ${calling_line}) couldn't open infile $infile\n";
        open(my $out, ">", $outname) or die "Subroutine $subname (called by script '${calling_script}', line ${calling_line}) couldn't create outfile $outname\n";

	# Loop over lines in the infile
	while(my $line = <$in>)
		{
		$line = text::trim($line);
		my @arr = split(",", $line);
		foreach my $el (@arr)	{	$el = text::trim($el);	}
		my $outline = join(",", @arr);
		print($out "$outline\n");
		}
	
        close($in);
        close($out);
	unlink($infile);
	move($outname, $infile)	
	}

return(1);

# end functions
