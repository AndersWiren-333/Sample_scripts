package fileTools;

# sub read_table($infile, $filetype)
# sub write_table($matrix_ref, $filetype, $outname)
# sub convert_table($infile, $infile_type, $outfile_type, $outfile_linebreaks_lin_win)
# sub write_list($array_ref, $filetype, $outname, $linebreaks)
# sub delete_table_columns($filename, $filetype, $col_index1, $col_index2 etc.)
# sub concat_PE($file1, $file2, $pattern)
# sub write_table_as_list($matrixref, $outname)
# sub trim_trinity_ID_patman($infile.pat, $outname)
# sub remove_multimatch_from_trinity_patman($infile.pat, $outname)
# sub read_diff_expr_sample_list($group_file)
# end sub list

########################################################## Universal perl module header ##########################################################

# perl_module_update

# Load libraries that this module depends on
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
my $thisfile = (__FILE__);
my $modules = "";
if($thisfile =~ m/^(.+)\//)	{	$modules = $1;	}
my $scripts = $modules;
$scripts =~ s/modules/scripts/;
my $maintain = $scripts;
$maintain =~ s/scripts/maintainance/;

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
require "$modules/combinatorics.pm";
require "$modules/db.pm";

# Create a timestamp string (can be attached to the name of logfiles, for example
my $timestamp = envir::timestamp();
my $rscript = "Rscript";

# end header

########################################################## Functions ##########################################################

# Reads a table from a file into a matrix. The file can be either a csv (comma separated values), a tsv (tab separated values) or a semicsv (semi-colon separated values) file
sub read_table
	{
	# Set error messages and accept input parameters
	my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
	my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$infile, \$filetype)'\n\nwhere".
	"\t\$infile is the file to be read\n".
	"\t\$filetype is the type of file to be read. Options are 'csv', 'tsv' and 'semicsv'\n\n";

	my @pars = @_ or die $usage;
	foreach my $el (@pars)  {       $el = text::trim($el);  }
	my $infile = shift @pars or die $usage;	# NB! This will not work if the argument is the number 0 (because it will then be interpreted as false). In that case you need to use 'shift(@pars) or 0', but it may lead to problems....
	my $filetype = shift @pars or die $usage;

	# Set field separator
	my $separator="";
	if($filetype eq "csv")	{	$separator = ",";	}
	elsif($filetype eq "tsv")  {       $separator = "\t";       }
	elsif($filetype eq "semicsv")  {       $separator = ";";       }
	else	{	die "\tUnknown filetype for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Options are 'csv, 'tsv' and 'semicsv'\n";	}

	# Open infile and declare empty matrix to hold results
	open(my $in, "<", $infile) or die "Subroutine $subname (called by script '${calling_script}', line ${calling_line}) couldn't open infile $infile\n";
	my @matrix=();

	while(my $line = <$in>)
		{
		$line = text::trim($line);
		if($line =~ /^\#/)	{	next;	}	# Skip any comment lines
		my @arr = split("$separator", $line);
		my @new_arr=();
		foreach my $el (@arr)
			{
			unless($el eq "")	{	$el = text::trim($el);	}
			$el =~ s/,/_/g;
			push(@new_arr, $el);
			}
		push(@matrix, [(@new_arr)]);
		}
	close($in);
	return(@matrix);
	}


# Prints a matrix (= table) to a file. The file can be either a csv (comma separated values), a tsv (tab separated values) or a semicsv (semi-colon separated values) file.
# If desired, several matrices can be written to the same file. References to additional matrices (except the first one) are given after the other arguments at the command line.
sub write_table
	{
	# Set error messages and accept input parameters
	my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
	my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$matrix_ref, \$filetype, \$outfile_name, \$linebreaks_lin_win, \$matrix_ref2, \$matrix_ref3 ...)'\n\nwhere".
	"\t\$matrix_ref is a reference to a matrix holding the data to be printed\n".
	"\t\$filetype is the type of file to be created. Options are 'csv', 'tsv' and 'semicsv'\n".
	"\t\$outfile_name is the name you want the outfile to have\n".
	"\t\$linebreaks_lin_win are the type of linebreaks to be used in the file. Options are 'lin' (Linux/Unix, \\n) and 'win' (Windows, \\r\\n)\n".
	"\t\$matrix_ref2, \$matrix_ref3 etc are references to additional matrices to be written to the same outfile, if desired.\n\n";
	my @pars = @_ or die $usage;
	foreach my $el (@pars)  {       $el = text::trim($el);  }
	my $matrixref = shift @pars or die $usage;      # NB! This will not work if the argument is the number 0 (because it will then be interpreted as false). In that case you need to use 'shift(@pars) or 0', but it may lead to proble$
	my $filetype = shift @pars or die $usage;
	my $outname = shift @pars or die $usage;
	my $linebreak = shift @pars or die $usage;
	my @all_matrices=();
	push(@all_matrices, $matrixref);
	if(@pars)	{	push(@all_matrices, @pars);	}

	# Set field separator
	my $separator="";
	if($filetype eq "csv")  {       $separator = ",";       }
	elsif($filetype eq "tsv")  {       $separator = "\t";       }
	elsif($filetype eq "semicsv")  {       $separator = ";";       }
	else    {       die "\tUnknown filetype for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Options are 'csv, 'tsv' and 'semicsv'\n";     }

	# Create outfile
	open(my $out, ">", $outname) or die "Subroutine $subname (called by script '${calling_script}', line ${calling_line}) couldn't create outfile $outname\n";

	# Loop over input matrices
	for(my $zz=0; $zz<=$#all_matrices; $zz++)
		{
		my $matref = $all_matrices[$zz];
		my @matrix = @{$matref};
	
		# Loop over rows in individual matrix
		for(my $cc=0; $cc<=$#matrix; $cc++)
			{
			my @arr = @{$matrix[$cc]};
			my $outline = join("$separator", @arr);
			if($linebreak eq "lin")	{	print($out "$outline\n");	}
			elsif($linebreak eq "win") {       print($out "$outline\r\n");       }
			}
			
		unless(scalar(@all_matrices)==1)	{	print $out "\n\n";	}
		}
	close($out);
	}

# Converts a table file (csv, tsv or semicsv) to another format (csv, tsv or semicsv)
sub convert_table
	{
	# Set error messages and accept input parameters
	my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
	my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$infile, \$infile_type, \$outfile_type, \$outfile_linebreaks_lin_win)'\n\nwhere".
	"\t\$infile is the file to be converted\n".
	"\t\$infile_type is the type of infile. Options are 'csv' (comma separated values), 'tsv' (tab separated values) and 'semicsv' (semi-colon separated values)\n".
	"\t\$outfile_type is the type of file to be created. Options are 'csv' (comma separated values), 'tsv' (tab separated values) and 'semicsv' (semi-colon separated values)\n".
	"\t\$linebreaks_lin_win are the type of linebreaks to be used in the created file. Options are 'lin' (Linux/Unix, \\n) and 'win' (Windows, \\r\\n)\n\n";

	my @pars = @_ or die $usage;
	foreach my $el (@pars)  {       $el = text::trim($el);  }

	my $infile = shift @pars or die $usage;      # NB! This will not work if the argument is the number 0 (because it will then be interpreted as false). In that case you need to use 'shift(@pars) or 0', but it may lead to proble$
	my $infile_type = shift @pars or die $usage;
	my $outfile_type = shift @pars or die $usage;
	my $linebreak = shift @pars or die $usage;

	my ($dir, $file, $basename, $suffix)=text::parse_path($infile);
	my $outname = "${basename}.${outfile_type}";

	# Read the infile
	my @matrix=fileTools::read_table($infile, $infile_type);

	# Create the outfile
	fileTools::write_table(\@matrix, $outfile_type, $outname, $linebreak);
	}

# Writes a list (an array) to a file
sub write_list
        {
        # Set error messages and accept input parameters
        my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
        my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$array_reference, \$outfile_name)'\n\nwhere".
        "\t\$array_reference is a refernce to the array that should be written\n".
	"\t\$outfile_name is the desired name of the outfile\n\n";

        my @pars = @_ or die $usage;
        foreach my $el (@pars)  {       $el = text::trim($el);  }

        my $arref = shift @pars or die $usage;      # NB! This will not work if the argument is the number 0 (because it will then be interpreted as false). In that case you need to use 'shift(@pars) or 0', but it may lead to proble$
        my $outname = shift @pars or die $usage;
	my @arr = @{$arref};

	open(my $out, ">", $outname) or die "Subroutine $subname (called by script '${calling_script}', line ${calling_line}) couldn't create outfile $outname\n";

	for(my $cc=0; $cc<=$#arr; $cc++)
		{
		print($out "$arr[$cc]\n");
		}
	
	close($out);
        }

# Deletes a specified column (or columns) from a table file (csv, tsv or semicsv). Retains both the original file and the new reduced file
# NB! The first column in a table has to be specified as "n" rather than "0", since Perl can't accept "0" as a command line argument.
sub delete_table_columns
        {
        # Set error messages and accept input parameters
        my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
        my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$filename, \$filetype, \$col_index1, \$col_index2 etc.)'\n\nwhere".
        "\t\$filename is the name of the file to delete columns from\n".
	"\t\$filetype is the type of file to be modified. Options are 'csv' (comma separated values), 'tsv' (tab separated values) and 'semicsv' (semi-colon separated values)\n".
        "\t\$col_index1 is the index of the first column to be deleted\n".
        "\t\$col_index2 is the index of the second column to be deleted, etc.\n".
        "\n".
        "\tNB! The first column in a table has to be specified as 'n' rather than '0', since Perl can't handle '0' as a command line argument.\n\n";
        my @pars = @_ or die $usage;
        foreach my $el (@pars)  {       $el = text::trim($el);  }
        my $infile = shift @pars or die $usage;      # NB! This will not work if the argument is the number 0 (because it will then be interpreted as false). In that case you need to use 'shift(@pars) or 0', but it may lead to proble$
	my $filetype = shift @pars or die $usage;
        my @cols = @pars or die $usage;

        # If @exclude_col_indices contains "n" (the represenation of zero), replace it with 0.
        foreach my $el (@cols)
                {
                if($el eq "n")  {       $el=0;       }
                }

	my $col_string = join("_", @cols);

        # Read the file into a matrix
	my @matrix=fileTools::read_table($infile, $filetype);

	# Delete selected columns
	my @new_matrix=matrixTools::delete_matrix_columns(\@matrix, @cols);

	# Print out the new file
	my $outname = "minus_${col_string}_${infile}";
	fileTools::write_table(\@new_matrix, $filetype, $outname, "lin");

        }

# Takes a list of paired end files (of some format, fastq, fasta, patman...) and concatenates them pairwise. The new, concatenated files are named based on excluding a specific text pattern that
# distinguishes file1 from file2 in a pair, from the name of file1. Returns an array with the new filenames. If the input or output files or both should be compressed after use, set parameter
# $zip to 'in', 'out' or 'both'.
sub concat_PE
	{
        # Set error messages and accept input parameters
        my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
        my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$namefile, \$pattern, \$zip, \$remove_infiles_y_n)'\n\nwhere".
        "\t\$namefile is a list of the names of the files that should be pairwise concatenated, pairs must be consecutive\n".
        "\t\$pattern is the text pattern that distinguishes the name of file1 from the name of file2 in a pair, e.g. '_R1_' or 'left'\n".
	"\t\$zip is an indicator of whether the input or output files or both should be compressed after use. Options are 'in', 'out' and 'both'.'\n".
	"\t\$remove_infiles_y_n is an indicator of whether the infiles should be deleted after concatenation. Options are 'y' and 'n'\n\n";

        my @pars = @_ or die $usage;
        foreach my $el (@pars)  {       $el = text::trim($el);  }
        my $namefile = shift @pars or die $usage;      # NB! This will not work if the argument is the number 0 (because it will then be interpreted as false). In that case you need to use 'shift(@pars) or 0', but it may lead to proble$
	my $pattern = shift @pars or die $usage;
	my $zip = shift @pars or die $usage;
	my $remove_infiles = shift @pars or die $usage;

	my @filenames=misc::read_list($namefile);
	my $num_files = scalar(@filenames);
	my $num_pairs = $num_files/2;

	my @new_names=();

	for(my $c=0; $c<$num_pairs; $c++)
        	{
        	my $infile1=shift(@filenames);
        	my $infile2=shift(@filenames);
	
		my ($file1_gz, $file1)=misc::check_compressed($infile1);
		my ($file2_gz, $file2)=misc::check_compressed($infile2);

		my $newname = $file1;
		$newname =~ s/$pattern/_PE-cat_/;
		system("cat $file1 $file2 > $newname");
		push(@new_names, $newname);

		if($remove_infiles eq "y")	{	unlink($file1, $file2);	}
		
		if(($remove_infiles eq "n") and ($zip eq "in"))
			{	
			gzip "$file1" => "${file1}.gz"; unlink($file1);
			gzip "$file2" => "${file2}.gz"; unlink($file2);
			}
		elsif(($remove_infiles eq "n") and ($zip eq "both"))
			{
			gzip "$file1" => "${file1}.gz"; unlink($file1);
			gzip "$file2" => "${file2}.gz"; unlink($file2);
			gzip "$newname" => "${newname}.gz"; unlink($newname);
			}
		if($zip eq "out")
			{	
			gzip "$newname" => "${newname}.gz"; unlink($newname);
			}
		}	
        return(@new_names);
	}

# Writes a table (matrix) to a file as a list (one element on each line, e.g. 'row1_ele1, $row1_ele2, $row2_ele1, $row2_ele2, $row3_ele1...')
sub write_table_as_list
	{
        # Set error messages and accept input parameters
        my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
        my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$matrix_reference, \$outname)'\n\nwhere".
        "\t\$matrix_reference is a reference to the matrix (table) to be written\n".
        "\t\$outname is the name that should be given to the resulting file\n\n";
        my @pars = @_ or die $usage;
        foreach my $el (@pars)  {       $el = text::trim($el);  }
        my $matref = shift @pars or die $usage;      # NB! This will not work if the argument is the number 0 (because it will then be interpreted as false). In that case you need to use 'shift(@pars) or 0', but it may lead to proble$
	my $outname = shift @pars or die $usage;
	my @matrix = @{$matref};
        open(my $out, ">", $outname) or die "Subroutine $subname (called by script '${calling_script}', line ${calling_line}) couldn't create outfile $outname\n";

	# Loop over lines in @matrix
	for(my $cc=0; $cc<=$#matrix; $cc++)
		{
		# Loop over elements in line and print each to the outfile
		my @arr = @{$matrix[$cc]};
		for(my $dd=0; $dd<=$#arr; $dd++)	{	print($out "$arr[$dd]\n");	}
		}

        close($out);
	}

# Trims fasta sequence ID lines produced by Trinity (transcriptome assembly software) in patman alignment output files
# Original format: "TRINITY_DN9_c0_g1_i1 len=268 path=[239:0-267] [-1, 239, -2]"
# Trimmed format: "TRINITY_DN9_c0_g1_i1"
sub trim_trinity_ID_patman
	{
        # Set error messages and accept input parameters
        my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
        my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$infile.pat, \$outname)'\n\nwhere".
        "\t\$infile.pat is the file (in patman output format = tsv) that should be processed\n".
	"\t\$outname is the desired name of the outfile\n\n";
        my @pars = @_ or die $usage;
        foreach my $el (@pars)  {       $el = text::trim($el);  }
        my $infile = shift @pars or die $usage;      # NB! This will not work if the argument is the number 0 (because it will then be interpreted as false). In that case you need to use 'shift(@pars) or 0', but it may lead to proble$
	my $outname = shift @pars or die $usage;

	# Open the infile
	open(my $in, "<", $infile) or die "Subroutine $subname (called by script '${calling_script}', line ${calling_line}) couldn't open infile $infile\n";
	open(my $out, ">", $outname) or die "Subroutine $subname (called by script '${calling_script}', line ${calling_line}) couldn't create outfile $outname\n";

	# Process each line one by one (slower than reading the whole file into memory and then processing(?) But conserves memory)
	while(my $line = <$in>)
		{
		$line = text::trim($line);
		my @arr = split("\t", $line);
		foreach my $el (@arr)	{	$el = text::trim($el);	}

		# Trim the ID
		if($arr[0] =~ /(TRINITY.+)( len=)/)	{	$arr[0] = $1;	}

		# Print out the new line
		my $outline = join("\t", @arr);
		print($out "$outline\n");
		}

        close($in);
        close($out);
	}
 
# If a read matches to more than one gene isoform in a transcriptome assembled with Trinity (which has a specific fasta header format that includes
# both gene- and isoform IDs), keep only the hit to the first isoform. This is useful if you later want to add the readcounts for each isoform into
# a total read count for the corresponding gene (if not, you may overestimate the readcount since the same physical read may count once for each isoform,
# when it must reasonably have come from only one of them)
sub remove_multimatch_from_trinity_patman
	{
        # Set error messages and accept input parameters
        my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
        my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$infile.pat \$outname)'\n\nwhere".
        "\t\$infile is the file to be processed\n".
        "\t\$outname is the desired name of the resulting file\n\n";
        my @pars = @_ or die $usage;
        foreach my $el (@pars)  {       $el = text::trim($el);  }
        my $infile = shift @pars or die $usage;      # NB! This will not work if the argument is the number 0 (because it will then be interpreted as false). In that case you need to use 'shift(@pars) or 0', but it may lead to proble$
	my $outname = shift @pars or die $usage;

	# Read the infile into a matrix
	my @matrix=fileTools::read_table($infile, "tsv");	

	# Loop over matrix and trim the isoform IDs away from the transcript IDs
	for(my $cc=0; $cc<=$#matrix; $cc++)
		{
		my $gene = $matrix[$cc][0];
		$gene =~ s/_i\d+//;
		$matrix[$cc][0] = $gene;
		}

	# Sort the matrix on transcript ID and read sequence
        @matrix = sort { $a->[0] cmp $b->[0] || $a->[1] cmp $b->[1] } @matrix;

	# Loop over matrix again to remove reads matching to multiple isoforms
	my @new_matrix=();
	my @comp_arr = @{$matrix[0]};
	my $comp_gene = $comp_arr[0];
	my $comp_read = $comp_arr[1];

	for(my $dd=1; $dd<=$#matrix; $dd++)
		{
		my @arr=@{$matrix[$dd]};
		my $gene = $arr[0];
		my $read = $arr[1];

		# If we are still on the same gene
		if($gene eq $comp_gene)
			{
			# If we are still on the same read
			if($read eq $comp_read)
				{	
				# If this is the last line of the matrix
				if($dd==$#matrix)	{	push(@new_matrix, [@comp_arr]);	}
				
				# If it is any other line
				else	{	next;	}
				}

			# If we are on a new read
			else
				{
				push(@new_matrix, [@comp_arr]);
				@comp_arr = @arr;
				$comp_read = $read;
				$comp_gene = $gene;

				# If this is the last line of the matrix
				if($dd== $#matrix)	{	push(@new_matrix, [@comp_arr]);	}
				}
			}

		# If we are on a new gene
		else
			{
			push(@new_matrix, [@comp_arr]);
			@comp_arr = @arr;
			$comp_read = $read;
			$comp_gene = $gene;

			# If this is the last line of the matrix
			if($dd==$#matrix)	{	push(@new_matrix, [@comp_arr]);	}
			}
		} # Loop over @matrix ends


	# Print the new matrix to an outfile
	fileTools::write_table(\@new_matrix, "tsv", $outname, "lin");

	# Count total and unique reads in the new file
	my @new_stats=misc::count_aligned_reads($outname, "nr");	

        return(@new_stats);
	}


# Reads a csv file with sample information needed when performing differential expression analysis on a gene expression matrix. The file lists
# the names of all samples that should be included in the analysis (not necessarily all samples in the dataset), along with information on
# sample groupings and which groups should be compared to each other for differential expression. The format should be this:
#
#	sample1,groupA
#	sample2,groupA
#	sample3,groupB
#	sample4,groupB
#	sample5,groupC
#	sample6,groupC
#	groupA_vs_groupB
#	groupA_vs_groupC
#	groupB_vs_groupC
#
# Returns two matrix references, one for a matrix holding the "sample-group" information and the other holding the "group_vs_group" information.
sub read_diff_expr_sample_list
	{
        # Set error messages and accept input parameters
        my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
        my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$group_file'\n\nwhere".
        "\t\$group_file is a csv file with sample information needed when performing differential expression analysis on a gene expression matrix.\n".
	"\tThe file lists the names of all samples that should be included in the analysis (not necessarily all samples in the dataset), along with information on\n".
	"\tsample groupings and which groups should be compared to each other for differential expression. The format should be this:\n\n".
	"\tsample1,groupA\n".
	"\tsample2,groupA\n".
	"\tsample3,groupB\n".
	"\tsample4,groupB\n".
	"\tsample5,groupC\n".
	"\tsample6,groupC\n".
	"\tgroupA_vs_groupB\n".
	"\tgroupA_vs_groupC\n".
	"\tgroupB_vs_groupC\n\n";

        my @pars = @_ or die $usage;
        foreach my $el (@pars)  {       $el = text::trim($el);  }
        my $infile = shift @pars or die $usage;      # NB! This will not work if the argument is the number 0 (because it will then be interpreted as false). In that case you need to use 'shift(@pars) or 0', but it may lead to proble$

	# Read the file into an aray
	my @lines=misc::read_list($infile);

	# Partition the contenst into two matrices
	my @sample_matrix=();
	my @comparisons=();

	for(my $cc=0; $cc<=$#lines; $cc++)
		{
		if($lines[$cc] =~ /_vs_/)
			{
			my @parts=split("_vs_", $lines[$cc]);
			push(@comparisons, [(@parts)])
			}
		else
			{
			my ($sample, $group) = split(",", $lines[$cc]);
			push(@sample_matrix, [($sample, $group)]);
			}
		}

	# Transform the @sample_matrix from 'sample1,groupA sample2,groupA' to 'groupA,sample1,sample2' format
	@sample_matrix = sort { $a->[1] cmp $b->[1] } @sample_matrix;
	my @group_matrix=();
	my $comp_group=333;	# Initializing with a number gives an error message if the value of the variable fails to be updated (to a string) later (=good, else we may never know it failed)
	my @templine=();

	for(my $dd=0; $dd<=$#sample_matrix; $dd++)
		{
		# If this is the first line in @sample_matrix
		if($dd==0)
			{
			$comp_group = $sample_matrix[$dd][1];
			push(@templine, $comp_group, $sample_matrix[$dd][0]);
			if($dd==$#sample_matrix)	{	die "There seems to be only one sample in your sample names file, so there is nothing to compare that sample to.\nAdd more samples and try again!\n";	}
			}

		else
			{
			# If we are still in the same group
			if($sample_matrix[$dd][1] eq $comp_group)
				{
				push(@templine, $sample_matrix[$dd][0]);
				
				# If this is the last line of @sample_matrix
				if($dd==$#sample_matrix)	{	push(@group_matrix, [(@templine)]);	}
				}

			# If we are in a new group
			else
				{
				push(@group_matrix, [@templine]);
				@templine=();
				$comp_group = $sample_matrix[$dd][1];
				push(@templine, $comp_group, $sample_matrix[$dd][0]);

				# If this is the last line of @sample_matrix
				if($dd==$#sample_matrix)        {       push(@group_matrix, [(@templine)]);       }
				}
			}
		} # Loop over @sample_matrix ends

	my @result=(\@group_matrix, \@comparisons);
	return(@result);
	}

# Compresses a file into a .gz file, in a platform independent way.
sub compress
	{
	# Set error messages
	my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
	my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$infile)'\n\nwhere".
	"\t\$infile is the file to be compressed\n";
	
	# Accept input parameters
	my @pars = @_ or die $usage;
	foreach my $el (@pars)  {       $el = text::trim($el);  }
	my $infile = shift @pars or die $usage;	
	my $newname = "${infile}.gz";
	
	gzip $infile => $newname;
	unlink($infile);
	}	
	
	
return(1);

# end functions
