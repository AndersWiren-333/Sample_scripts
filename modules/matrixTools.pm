package matrixTools;

# sub testprint_matrix(\@matrix)
# sub transpose_matrix(\@matrix)
# sub generate_rand_matrix($rows, $columns, $min, $max, $decimals)
# sub create_zero_matrix($num_rows, $num_cols)
# sub merge_matrices_col(\@matrix1)
# sub split_matrix_on_col_value(\@matrix, $col_index)
# sub get_matrix_columns_as_rows($matrixref, $col_index1, $col_index2...)
# sub delete_matrix_columns($matrixref, $col_index1, $col_index2...)
# sub add_length_to_matrix_gene_id($matrixref, $length_col_index, $header_y_n)
# sub add_cols_to_matrix($matrixref, $col1_ref, $col2_ref etc...)
# sub sync_genelists_on_sortkey($matrix1.csv, $sortkey1, $matrix2.csv, $sortkey2)
# sub apply($matrixref, $row_or_col, $function)
# sub dim($matrix_reference)
# sub replace_matrix_cols(\@matrix, \@replacement_cols_matrix, $col_to_replace1, $col_to_replace2 etc...)
# sub get_matrix_columns($matrix_reference, $col1, $col2 etc...)
# sub matrix_equality($matrixref1, $matrixref2)
# sub add_seq_lengths_to_matrix_from_fasta($file_or_var, $filename_or_matrixref, $header_y_n, $fasta_file)
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

# Prints a given matrix on screen
sub testprint_matrix
	{
	# Set error messages
	my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
	my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$matrixref)'\n\nwhere".
	"\t\$matrixref is a reference to the matrix to be printed\n\n";
	
	# Accept input parameters
	my @pars = @_ or die $usage;
	foreach my $el (@pars)  {       $el = text::trim($el);  }
	my $matrixref = shift @pars or die $usage;	# NB! This will not work if the argument is the number 0 (because it will then be interpreted as false). In that case you need to use 'shift(@pars) or 0', but it may lead to problems....
	
	# Process matrix
	my @matrix = @{$matrixref};

	# Loop over rows in matrix
	for(my $c=0; $c<=$#matrix; $c++)
		{
		my @arr=@{$matrix[$c]};
				
		# Loop over columns in row
		for(my $d=0; $d<=$#arr; $d++)	{	print("$arr[$d]\t");	}
		print("\n");
		}
	}

# Transposes a given matrix (converts rows to columns and columns to rows)
sub transpose_matrix
	{
	my $usage = "Syntax error for sub transpose_matrix. Correct usage: 'matrixTools::transpose_matrix(\@matrix)'\n";
	my $matrix_ref = $_[0] or die $usage;
	my @matrix = @{$matrix_ref};
	my $num_rows=scalar(@matrix);
	my @first_row=@{$matrix[0]};
	my $num_cols=scalar(@first_row);

	my @t_matrix=();

	# Loop over columns in original matrix
	for(my $cc=0; $cc<$num_cols; $cc++)
		{
		# Loop over rows in riginal matrix
		for(my $dd=0; $dd<$num_rows; $dd++)
			{
			$t_matrix[$cc][$dd] = $matrix[$dd][$cc];
			}
		}
	return(@t_matrix);
	}

# Generates a matrix with given dimensions, filled with random numbers between given values and with a given number of decimals
sub generate_rand_matrix
	{
	my $usage = "Syntax error for sub generate_rand_matrix. Correct usage: 'matrixTools::generate_rand_matrix(\$rows, \$columns, \$min, \$max, \$decimals);'\n";
	my $rows = $_[0] or die $usage;
	my $cols = $_[1] or die $usage;
	my $min = $_[2] or die $usage;
	my $max = $_[3] or die $usage;
	my $dec = $_[4] or 0;
	my @matrix=();

	for(my $cc=0; $cc<$rows; $cc++)
		{
		for(my $dd=0; $dd<$cols; $dd++)
			{
			$matrix[$cc][$dd] = stats::rand_between($min, $max, $dec);
			}
		}
	return(@matrix);
	}

# Creates a matrix with the specified numbers of rows and columns, filled with zeros
sub create_zero_matrix
	{
	my $usage="\nSyntax error for sub create_zero_matrix. Correct usage: 'matrixTools::create_zero_matrix(\$num_rows, \$num_cols)'\n\nwhere".
        "\t\$num_rows is the desired number of rows\n".
        "\t\$num_cols is the desired number of columns\n";
	my @pars = @_ or die $usage;
	foreach my $el (@pars)	{	$el = text::trim($el);	}
	my $nrows = $pars[0] or die $usage;
	my $ncols = $pars[1] or die $usage;
	my @matrix=();

	# Loop over rows
	for(my $cc=0; $cc<$nrows; $cc++)
		{
		# Loop over columns
		for(my $dd=0; $dd<$ncols; $dd++)
			{
			$matrix[$cc][$dd]=0;
			}		
		}

	return(@matrix);
	}

# Merges two or more matrices (tables) in csv format into one. The function checks whether the same entry in the first column of one matrix has an
# identical counterpart in the first column of the other matrices (this may for example be a gene id), and if so places all columns from all matrices
# corresponding to that entry on the same row in the output matrix. Entries that don't have a value in a particular matrix will be assigned a value of "0"
# for the respective column.
sub merge_matrices_col
	{
	my $usage="\nSyntax error for sub merge_matrices_col. Correct usage: 'matrixTools::merge_matrices_col(\$header_y_n, \$outfile_name.csv, \$matrix_1.csv, \$matrix_2.csv, \$matrix_3.csv etc.)'\n\nwhere".
	"\t\$header_y_n is an indicator of whether the infiles have headers or not. Options are 'y' and 'n'\n".
	"\t\$outfile_name.csv is the desired name of the outfile\n".
        "\t\$matrix_1.csv, \$matrix_2.csv etc. are the matrix files (in csv format) to be merged\n\n";
	my @pars = @_ or die $usage;
	foreach my $el (@pars)	{	$el = text::trim($el);	}
	my $headers_y_n = shift(@pars);
	my $outname = shift @pars;
	my @infiles = @pars;
	my @outfile_headers=();
	my @ids=();
	my @matrices=();

	# Loop over infiles
	for(my $cc=0; $cc<=$#infiles; $cc++)
		{
		# Read in the matrix and format headers for the outfile
		my $matrix_name = $infiles[$cc];
		my @matrix_all = misc::file_to_matrix($matrix_name, "comma", "all", "y");
		my @headers = @{$matrix_all[1]};
		my @matrix = @{$matrix_all[0]};
		@matrix = sort { $a->[0] cmp $b->[0] } @matrix;
		if($cc==0)	{	push(@outfile_headers, @headers);	}
		else	{	shift(@headers); push(@outfile_headers, @headers);	}

		# Loop over matrix and add all gene ids to the list of ids
		for(my $dd=0; $dd<=$#matrix; $dd++)
			{
			my $id = $matrix[$dd][0];
			push(@ids, $id);
			}

		# Add a reference to the matrix to @matrices
		push(@matrices, \@matrix);
		}

	my @uni_ids=misc::unique_list(\@ids, "alph");
	my @grand_matrix = matrixTools::create_zero_matrix(scalar(@uni_ids), scalar(@outfile_headers));
	my $col_index=1;

	# Set the unique feature ID's as first column of @grand_matrix
	for(my $ee=0; $ee<=$#grand_matrix; $ee++)	{	$grand_matrix[$ee][0] = $uni_ids[$ee];	}

	# For each matrix, check each ID in @uni_ids and fill in its expression value in the @grand_matrix
	for(my $ff=0; $ff<=$#matrices; $ff++)
		{
		my @matrix = @{$matrices[$ff]};
		@matrix = sort { $a->[0] cmp $b->[0] } @matrix;
		my $matrix_cols=scalar(@{$matrix[0]});	
		my $ref_ind=0;
		my $mat_ind=0;

		# Loop over @grand_matrix
		while($ref_ind<=$#grand_matrix)
			{
			# If reference id is before matrix id, go to next reference id
			if($grand_matrix[$ref_ind][0] lt $matrix[$mat_ind][0])	{	$ref_ind++;	}
	
			# Else, if reference id and matrix id are the same, assign expression values to @grand_matrix
			elsif($grand_matrix[$ref_ind][0] eq $matrix[$mat_ind][0])
				{
				# Add each expression value for the gene in question in the matrix in question to @grand_matrix
				for(my $col=1; $col<$matrix_cols; $col++)
					{
					my $grand_col = $col_index+$col-1;
					$grand_matrix[$ref_ind][$grand_col] = $matrix[$mat_ind][$col];
					}
				$ref_ind++;
				$mat_ind++;
				}
			# Else, if reference id is after matrix id, go to next reference id
			elsif($grand_matrix[$ref_ind][0] gt $matrix[$mat_ind][0])	{	$ref_ind++;	}
			if($mat_ind > $#matrix)	{	$mat_ind--;	}
			}
		$col_index = $col_index+$matrix_cols-1;
		}

	unshift(@grand_matrix, \@outfile_headers);
	matrixTools::matrix_to_csv(\@grand_matrix, $outname);	
	return(@grand_matrix);
	}


# Splits a matrix on rows to create a set of new matrices each representing the rows in the original matrix that have a unique value of the given column.
# For example @submatrix_alma has all the rows from the original matrix that have the value "alma" in column number 3, and @submatrix_kalle those rows where
# column 3 has a value of "kalle". Except for specifying the matrix name and the column number, also specify if the column has numeric ("num") or alphabetic ("alph") values
sub split_matrix_on_col_value
	{
	# Set error messages and accept input parameters
        my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
        my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$matrixref, \$num_or_alph, \$col_index)'\n\nwhere".
        "\t\$matrixref is a reference to the matrix to be split\n".
	"\t\$num_or_alph is an indicator of whether the column of interest holds numeric or alphabetic values (options: 'num' or 'alph')\n".
	"\t\$col_index is the index of the column in the matrix that the splitting should be based on\n";

        my @pars = @_ or die $usage;
        foreach my $el (@pars)  {       $el = text::trim($el);  }
        my $matrixref = shift @pars or die $usage;      # NB! This will not work if the argument is the number 0 (because it will then be interpreted as false). In that case you need to use 'shift(@pars) or 0', but it may lead to proble$
	my $num_or_alph = shift @pars or die $usage;
	my $col = shift @pars or 0;
	my @matrix=@{$matrixref};

	# Put all the values of the chosen column into an array
	my @col_values=();
	for(my $cc=0; $cc<=$#matrix; $cc++)	{	push(@col_values, $matrix[$cc][$col]);	}

	# Sort matrix on chosen column
	if($num_or_alph eq "alph")	{	@matrix = sort { $a->[$col] cmp $b->[$col] } @matrix;	}
	elsif($num_or_alph eq "num")	{	@matrix = sort { $a->[$col] <=> $b->[$col] } @matrix;	}

	# Find out the number and identity of unique values in that column	
	my @unique_values=misc::unique_list(\@col_values, $num_or_alph);
	my $num_values=scalar(@unique_values);

	# Loop over rows in @matrix, and for each new value in column $col, start printing lines to a new sub-matrix
	my @matrix_set=();
	my $sub_matrix_index=0;
	my $sub_matrix_line_index=0;
	
	for(my $dd=0; $dd<=$#matrix; $dd++)
		{
		my $sub_matrix_value = $unique_values[$sub_matrix_index];
		my $value = $matrix[$dd][$col];
		
		if($value eq $sub_matrix_value)
			{
			$matrix_set[$sub_matrix_index][$sub_matrix_line_index] = $matrix[$dd]; 
			$sub_matrix_line_index++;
			}
		elsif($value ne $sub_matrix_value)
			{
			$sub_matrix_index++;
			$sub_matrix_line_index=0;

			$matrix_set[$sub_matrix_index][$sub_matrix_line_index] = $matrix[$dd];
			$sub_matrix_line_index++;
			}
		}
	
        return(@matrix_set);
	}

# Prints a given 3d matrix on the screen (given a reference to it)
sub testprint_3d_matrix
	{
        # Set error messages and accept input parameters
        my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
        my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$3d_matrixref)'\n\nwhere".
        "\t\$3d_matrixref is a reference to the 3d matrix to be printed\n\n";

        my @pars = @_ or die $usage;
        foreach my $el (@pars)  {       $el = text::trim($el);  }
        my $matrixref_3d = shift @pars or die $usage;      # NB! This will not work if the argument is the number 0 (because it will then be interpreted as false). In that case you need to use 'shift(@pars) or 0', but it may lead to proble$
	
	my @matrix_3d = @{$matrixref_3d};

        for(my $cc=0; $cc<=$#matrix_3d; $cc++)
                {
                print("This is matrix number ${cc}:\n\n");
                my @matrix_2d = @{$matrix_3d[$cc]};

                for(my $dd=0; $dd<=$#matrix_2d; $dd++)
                        {
                        my @arr = @{$matrix_2d[$dd]};

                        for(my $ee=0; $ee<=$#arr; $ee++)
                                {
                                print("$arr[$ee]\t");
                                }
                        print("\n");
                        }
                print("---------------------------------------------\n");
		}
	}

# Retrieves specific columns (one or many) from a matrix, specified by the index (number) of those columns in the matrix. If one column is requested,
# the function returns it as an array. If more columns are requested, they are returned as a matrix where each row is a reference to an array holding
# a specific column. The first column is called 1, the second 2 etc (i.e. not standard perl indexing).
sub get_matrix_columns_as_rows
	{
	# Set error messages and accept input parameters
	my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
	my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$matrixref, \$col_index1, \$col_index2 etc)'\n\nwhere".
	"\t\$matrixref is a refernce to the matrix from which columns should be retrieved\n".
	"\t\$col_index1 is the index of the first column to be retrieved. NB! The first column is called '1' rather than '0'.\n".
	"\t\$col_index2 is the index of the second column to be retrieved, etc.\n\n";
	my @pars = @_ or die $usage;
	foreach my $el (@pars)  {       $el = text::trim($el);  }
	my $matrixref = shift @pars or die $usage;
	my @cols = @pars or die $usage;
	@cols=misc::shift_input_cols(\@cols);

	# Processing
	my $num_cols = scalar(@cols);
	my @matrix = @{$matrixref};

	# Get the columns
	my @one=();
	my @all=();

	# Loop over rows in matrix
	for(my $row=0; $row<=$#matrix; $row++)
		{
		# Loop over the requested columns
		for(my $col=0; $col<=$#cols; $col++)
			{
			my $index = $cols[$col];
			my $value = $matrix[$row][$index];
			if($num_cols>1)	{	$all[$col][$row] = $value;	}
			elsif($num_cols==1) {       $one[$row] = $value;      }
			}
		}

	if($num_cols>1)	{	return(@all);	}
	elsif($num_cols==1) {       return(@one);   }
	}

# Deletes specified columns from a matrix. NB! The first column in a matrix is called 1 rather than 0 (i.e. not perl native indexing).
sub delete_matrix_columns
	{
	# Set error messages and accept input parameters
	my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
	my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$matrixref, \$col_index1, \$col_index2 etc)'\n\nwhere".
	"\t\$matrixref is a reference to the matrix from which columns should be deleted\n".
	"\t\$col_index1 is the index of the first column to be deleted. NB! The first column in a matrix is called 1 rather than 0.\n".
	"\t\$col_index2 is the index of the second column to be deleted, etc.\n\n";

	my @pars = @_ or die $usage;
	foreach my $el (@pars)  {       $el = text::trim($el);  }
	my $matrixref = shift @pars or die $usage;      # NB! This will not work if the argument is the number 0 (because it will then be interpreted as false). In that case you need to use 'shift(@pars) or 0', but it may lead to proble$
	my @exclude_col_indices = @pars or die $usage;
	my @matrix = @{$matrixref};
	
	@exclude_col_indices=misc::shift_input_cols(\@exclude_col_indices);

	my @new_matrix=();

	# Loop over rows in @matrix
	for(my $cc=0; $cc<=$#matrix; $cc++)
		{
		my @arr = @{$matrix[$cc]};	
		my $reduce_factor=0;

		# For each column that should be removed, splice it out from the array
		for(my $dd=0; $dd<=$#exclude_col_indices; $dd++)
			{
			my $index = $exclude_col_indices[$dd];
			$index -= $reduce_factor;
			splice(@arr, $index, 1);
			$reduce_factor++;
			}
		push(@new_matrix, \@arr);
		}

	return(@new_matrix);
	}

# Takes a gene expression with a gene length column and adds the length information to the ID/name column of every gene, e.g. "Etteplan_B" becomes "Etteplan_B_1_237"
sub add_length_to_matrix_gene_id
	{
	# Set error messages and accept input parameters
	my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
	my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$matrixref, \$length_col_index, \$header_y_n)'\n\nwhere".
	"\t\$matrixref is a reference to the matrix to be processed\n".
	"\t\$length_col_index is the column number (minus one) of the column holding the length information. If this is 0, set it to 'n' (Perl can't handle zero as a command line argument)\n".
	"\t\$header_y_n is an indicator of whether the matrix has a header or not. 'y' is yes and 'n' is no\n\n";
	my @pars = @_ or die $usage;
	foreach my $el (@pars)  {       $el = text::trim($el);  }
	my $matrixref = shift @pars or die $usage;      # NB! This will not work if the argument is the number 0 (because it will then be interpreted as false). In that case you need to use 'shift(@pars) or 0', but it may lead to proble$
	my $col_ind = shift @pars or die $usage;
	my $header_y_n = shift @pars or die $usage;

	if($col_ind eq "n")	{	$col_ind = 0;	}
	my @matrix = @{$matrixref};
	my $start_at="Faster_Tinne";
	if($header_y_n eq "y")	{	$start_at = 1;	}
	elsif($header_y_n eq "n")	{	$start_at = 0;	}

	# Loop over @matrix and add length information to the gene/feature id
	for(my $cc=$start_at; $cc<=$#matrix; $cc++)
		{
		my $length = $matrix[$cc][$col_ind];
		my $new_id = "$matrix[$cc][0]"."_1_"."$length";
		$matrix[$cc][0] = $new_id;
		}

	return(@matrix);
	}

# Adds one or more columns to a matrix. The new columns are specified as array references.
sub add_cols_to_matrix
	{
	# Set error messages and accept input parameters
	my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
	my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$matrixref, \$colref_1, \$colref_2 etc...)'\n\nwhere".
	"\t\$matrixref is a reference to the matrix that columns should be added to\n".
	"\t\$colref_1 is a reference to a column containing the values that should become the first new column\n".
	"\t\$colref_2 is a reference to a column containing the values that should become the secons new column\n\n";
	my @pars = @_ or die $usage;
	foreach my $el (@pars)  {       $el = text::trim($el);  }
	my $matrixref = shift @pars or die $usage;	# NB! This will not work if the argument is the number 0 (because it will then be interpreted as false). In that case you need to use 'shift(@pars) or 0', but it may lead to problems....
	my @columns = @pars or die $usage;

	my @matrix = @{$matrixref};
	
	# If the matrix is empty this function will not work, so we have to create a dummy row to start with
	my @first_row = (1,2,3);
	my @second_row = (4,5,6);
	my $empty="n";
	if(scalar(@matrix)==0)	{	$matrix[0] = (\@first_row);	$matrix[1] = (\@second_row); $empty="y";	}
	
	# Transpose matrix
	my @t_matrix=matrixTools::transpose_matrix(\@matrix);
	
	# Add new columns (as rows, since the matrix is now transposed)
	foreach my $el (@columns)	{	push(@t_matrix, $el);	}
	
	# If the input matrix was empty, remove the now three (before transposition two) dummy rows we added earlier to get something to add columns to
	if($empty eq "y")	{	shift(@t_matrix); shift(@t_matrix); shift(@t_matrix);	}
	
	# Transpose back
	my @new_matrix=matrixTools::transpose_matrix(\@t_matrix);
		
	return(@new_matrix);	
	}

# Takes two genelists and removes any genes that are only present in one of them, i.e. it retains those genes that are common to the two lists
# (i.e. it produces the intersection of the lists, but doesn’t join them – it outputs two new lists). For this to be possible, the lists need
# to use the same type of gene identifiers, and the column number of these identifiers are given as arguments at the command line.
sub sync_genelists_on_sortkey
	{
	# Set error messages
	my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
	my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$matrix1.csv, \$sortkey1, \$matrix2.csv, \$sortkey2)'\n\nwhere".
	"\t\$matrix1.csv is the first genelist\n".
	"\t\$sortkey1 is the column number in the first genelist where the common identifiers are stored\n".
	"\t\tNB! NB! If this should be 0, set it to 'zero' (perl can't handle the number 0 as a command line argument)\n".
	"\t\$matrix2.csv is the second genelist\n".
	"\t\$sortkey2 is the column number in the second genelist where the common identifiers are stored\n".
	"\t\tNB! NB! If this should be 0, set it to 'zero' (perl can't handle the number 0 as a command line argument)\n\n";

	# Accept input parameters
	my @pars = @_ or die $usage;
	foreach my $el (@pars)  {       $el = text::trim($el);  }
	
	my $infile_1 = shift(@pars) or die "Parameter error \$infile_1\n\n".$usage;
	my $sortkey_1 = shift(@pars) or die "Parameter error \$sortkey_1\n\n".$usage;
	my $infile_2 = shift(@pars) or die "Parameter error \$infile_2\n\n".$usage;
	my $sortkey_2 = shift(@pars) or die "Parameter error \$sortkey_2\n\n".$usage;

	if($sortkey_1 eq "zero")	{	$sortkey_1=0;	}
	if($sortkey_2 eq "zero")	{	$sortkey_2=0;	}
	
	my $outname_1 = "synced_${infile_1}";
	my $outname_2 = "synced_${infile_2}";

	# Read infiles
	my @matrix1=fileTools::read_table($infile_1, "csv");
	my @matrix2=fileTools::read_table($infile_2, "csv");

	# Remove headers
	my $headerref1 = shift(@matrix1);
	my $headerref2 = shift(@matrix2);
	
	my @new_matrix1=();
	my @new_matrix2=();
	
	# Sort input matrices for faster processing
	my @s_matrix1 = sort { $a->[$sortkey_1] cmp $b->[$sortkey_1] } @matrix1;
	my @s_matrix2 = sort { $a->[$sortkey_2] cmp $b->[$sortkey_2] } @matrix2;

	# Loop over genes in first matrix
	FIRST: for(my $c=0; $c<=$#s_matrix1; $c++)
		{
		# Loop over genes in second matrix
		SECOND: for(my $d=0; $d<=$#s_matrix2; $d++)
			{
			if($s_matrix1[$c][$sortkey_1] eq $s_matrix2[$d][$sortkey_2])
				{
				push(@new_matrix1, $s_matrix1[$c]);
				push(@new_matrix2, $s_matrix2[$d]);
				next FIRST;
				}
			}
		}

	# Add headers to new matrices
	unshift(@new_matrix1, $headerref1);
	unshift(@new_matrix2, $headerref2);
		
	# Write new matrices to files
	fileTools::write_table(\@new_matrix1, "csv", $outname_1, "lin");
	fileTools::write_table(\@new_matrix2, "csv", $outname_2, "lin");
	}

# Applies a function (mean, median, min, max) to all rows ($row_or_col=”row”) or columns ($row_or_col=”col”)
# in a matrix and returns the resulting values as an array.
sub apply
	{
	# Set error messages
	my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
	my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$matrixref, \$row_or_col, \$function)'\n\nwhere".
	"\t\$matrixref is a reference to a matrix\n".
	"\t\$row_or_col  Set to 'row' if '\$function' should be applied to rows, or to 'col' if it should be applied to columns\n".
	"\t\$function is the function that should be applied to the rows or columns. Options are 'mean', 'median', 'min', 'max', 'geomean', 'stdev', 'variance'\n\n";
	
	# Accept input parameters
	my @pars = @_ or die $usage;
	foreach my $el (@pars)  {	$el = text::trim($el);	}
	my $matrixref = shift @pars or die $usage;	# NB! This will not work if the argument is the number 0 (because it will then be interpreted as false). In that case you need to use 'shift(@pars) or 0', but it may lead to problems....
	my $row_or_col = shift @pars or die $usage;
	my $function = shift @pars or die $usage;

	my @matrix = @{$matrixref};
	my @final_results=();

	if($row_or_col eq "col")	{	@matrix=matrixTools::transpose_matrix(\@matrix);	}
	
	# Loop over rows (or columns) in matrix and compute aggregate values
	my @aggregate_values=();
	my $result="";
	for(my $cc=0; $cc<=$#matrix; $cc++)
		{
		if($function eq "mean")	{	$result=stats::mean($matrix[$cc]);	}
		elsif($function eq "median")	{	$result=stats::median($matrix[$cc]);	}
		elsif($function eq "min")	{	$result=stats::min($matrix[$cc]);	}
		elsif($function eq "max")	{	$result=stats::max($matrix[$cc]);	}
		elsif($function eq "geomean")	{	$result=stats::geometric_mean($matrix[$cc]);	}
		elsif($function eq "stdev")	{	$result=stats::stdev($matrix[$cc]);	}
		elsif($function eq "variance")	{	$result=stats::variance($matrix[$cc]);	}
		else	{	die "Function $function isn't supported by this subroutine ${subname}. $usage";	}
		
		push(@aggregate_values, $result);
		}
		
	return(@aggregate_values);
	}

# Returns the number of rows and columns of a matrix (given as a matrix reference), as an array. 
sub dim
	{
	# Set error messages
	my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
	my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$matrixref)'\n\nwhere".
	"\t\$matrixref is a reference to a matrix\n\n";
	
	# Accept input parameters
	my @pars = @_ or die $usage;
	foreach my $el (@pars)  {       $el = text::trim($el);  }
	my $matrixref = shift @pars or die $usage;	# NB! This will not work if the argument is the number 0 (because it will then be interpreted as false). In that case you need to use 'shift(@pars) or 0', but it may lead to problems....
	
	# Processing
	my @matrix = @{$matrixref};
	my $num_rows = scalar(@matrix);
	my @transposed_matrix=matrixTools::transpose_matrix(\@matrix);
	my $num_cols = scalar(@transposed_matrix);
	my @dimensions = ($num_rows, $num_cols);
	
	return(@dimensions);
	}
	
# Takes a matrix and a range of rows and columns and replaces that range (chunk) with another specified marix (of those dimensions).
# Specify the range to be replaced as “1_4_to_3_7” (replace the chuck starting at row 1, column 4 and ending at row 3, column 7).
sub replace_matrix_chunk
	{
	# Set error messages
	my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
	my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$argument1, \$argument2, \@arguments3)'\n\nwhere".
	"\t\$argument1 can be either 'X' or 'Y'\n".
	"\t\$argument2 is a file in XX format\n".
	"\t\@arguments3 is a list of arguments\n\n";
	
	# Accept input parameters
	my @pars = @_ or die $usage;
	foreach my $el (@pars)  {       $el = text::trim($el);  }
	my $matrixref = shift @pars or die $usage;
	my $dim_string = shift @pars or die $usage;
	my $chunkref = shift @pars or die $usage;

	# Processing
	my @matrix = @{$matrixref};
	my @chunk = @{$chunkref};
	
	my ($startrow, $startcol, $to, $endrow, $endcol) = split("_", $dim_string);
	$startrow--; $endrow--; $startcol--; $endcol--;
	
	# If the specified chunk of the matrix and the replacement chunk don't have equal dimensions,
	# terminate the program with an error message.
	my @specdims = (($endrow-$startrow+1),($endcol-$startcol+1));
	my @repdims=matrixTools::dim(\@chunk);
	my $equal_dims=listTools::array_equality(\@specdims, \@repdims);
	if($equal_dims eq "n")	{	die "The specified part of the matrix doesn't have the same size as the replacement chunk.\nCheck the dimensions and try again!\n";	}
	
	#print "$specdims[0] $specdims[1]\n";
	#print "$repdims[0] $repdims[1]\n";
	
	# Loop over specified row and columns in @matrix
	for(my $mrow=$startrow; $mrow<=$endrow; $mrow++)
		{
		my $rrow = $mrow - $startrow;
		for(my $mcol=$startcol; $mcol<=$endcol; $mcol++)
			{
			my $rcol = $mcol - $startcol;
			$matrix[$mrow][$mcol] = $chunk[$rrow][$rcol];
			}
		}
	
	return(@matrix);
	}

# Retrieves specific columns (one or many) from a matrix, specified by the index (number) of those columns in the matrix. They are returned as a matrix.
# Note that more columns than one must be requested (otherwise use subroutine ‘get_matrix_columns_as_rows’). Column numbering begins at 1 rather than 0.
sub get_matrix_columns
	{
	# Set error messages
	my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
	my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$matrix_reference, \$column_index1, \@column_index2 etc...)'\n\nwhere".
	"\t\$matrix_reference is a reference to the matrix from which columns should be retrieved\n".
	"\t\$column_index1, column_index2 etc ... are the numbers from left to right of the columns that should be retrieved (numbering begins at 1)\n\n";
	
	# Accept input parameters
	my @pars = @_ or die $usage;
	foreach my $el (@pars)  {       $el = text::trim($el);  }
	my $matref = shift @pars or die $usage;	
	my @cols = @pars or die $usage;
	
	# Processing
	my @matrix_rows=matrixTools::get_matrix_columns_as_rows($matref, @cols);
	my @transposed_matrix=matrixTools::transpose_matrix(\@matrix_rows);

	return(@transposed_matrix);
	}
	
	
# Checks whether two matrices (specified as matrix references) are equal, and returns “y” if they are and “n” if they aren’t.
sub matrix_equality
	{
	# Set error messages
	my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
	my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$matrix_reference1, \$matrix_reference2)'\n\nwhere".
	"\t\$matrix_reference1 is a reference to the first matrix to be compared\n".
	"\t\$matrix_reference2 is a reference to the second matrix to be compared\n\n";
	
	# Accept input parameters
	my @pars = @_ or die $usage;
	foreach my $el (@pars)  {       $el = text::trim($el);  }
	my $matref1 = shift @pars or die $usage;
	my $matref2 = shift @pars or die $usage;
	my @matrix1 = @{$matref1};
	my @matrix2 = @{$matref2};
	
	# Check if the matrices have the same dimensions (if they don't they are not equal)
	my @dim1=matrixTools::dim($matref1);
	my @dim2=matrixTools::dim($matref2);
	my $same_dims=listTools::array_equality(\@dim1, \@dim2);
	my $matrices_equal = "y";
	if($same_dims eq "n")	{	$matrices_equal = "n";	}
	
	# If the dimensions are equal, check if all values are equal
	elsif($same_dims eq "y")
		{
		# Loop over rows in matrix 1
		for(my $cc=0; $cc<=$#matrix1; $cc++)
			{
			my @row1 = @{$matrix1[$cc]};
			my @row2 = @{$matrix2[$cc]};
			my $rows_equal=listTools::array_equality(\@row1, \@row2);
			if($rows_equal eq "n")	{	$matrices_equal = "n"; last;	}
			my $rad = 1+$cc;
			print "Rows $rad equal: $rows_equal\n";
			}
		}

	else	{	die "\n\tEquality of matrices couldn't be evalueated because of an internal error\n";	}
	return($matrices_equal);
	}

# Takes a matrix (as a csv file – set $file_or_var to ‘file’ - or matrix reference – set $file_or_var to ‘var’ - where
# each row represents a genomic feature (e.g. a gene or transcript and a specified column is a feature ID)) and a fasta
# file containing the DNA/RNA/protein sequence of the features in the matrix. Using the fasta file the subroutine computes
# the length of each feature in the matrix, and attaches these lengths to the matrix as a new column. Set $header_y_n to ‘y’
# if the original matrix has a header, otherwise to ‘n’. If the matrix was read from a csv file, the new matrix will also
# be printed to a csv file (as well as being returned as a matrix).
sub add_seq_lengths_to_matrix_from_fasta
	{
	# Set error messages and accept input parameters
	my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
	my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$file_or_var, \$filename_or_matrixref, \$id_col, \$header_y_n, \$fasta_file)'\n\nwhere".
	"\t\$file_or_var should be set to 'var' if the input is a reference to a matrix, or 'file' if the input is a csv file (and the output also should be)\n".
	"\t\$filename_or_matrixref is either the name of the input csv file or a reference to the matrix to be processed\n".
	"\t\$id_col is a numberr representing the column in the matrix that holds feature IDs. Numbering starts at 1 (i.e. not native perl indexing)\n".
	"\t\$header_y_n - set this to 'y' if the input file/matrix has a header or 'n' if it doesn't\n".
	"\t\$fasta_file is the fasta file from which lengths should be computed\n\n";
	
	my @pars = @_ or die $usage;
	foreach my $el (@pars)  {       $el = text::trim($el);  }
	my $file_or_var = shift @pars or die $usage;
	my $matrix_file_or_ref = shift @pars or die $usage;
	my $id_col = shift @pars or die $usage;
	$id_col--;	# Converts this column index to perl native indexing (which starts at 0, which can't be used as a command line argument)
	my $header_y_n = shift @pars or die $usage;
	my $fasta = shift @pars or die $usage;
	
	my @matrix=();
	if($file_or_var eq "file")	{	@matrix=fileTools::read_table($matrix_file_or_ref, "csv");	}
	elsif($file_or_var eq "var")	{	@matrix = @{$matrix_file_or_ref};	}
	else	{	die "\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Parameter \$file_or_var must be either 'file' or 'var'\n\n";	}
	
	# If there is a header, remove it from the matrix and store it in an array
	my @header=();
	if($header_y_n eq "y")
		{
		my $headerref = shift(@matrix);
		@header = @{$headerref};
		push(@header, "length");
		}

	# Get sequence lengths from fasta file
	my @length_matrix=fastaTools::get_seq_lengths_from_fasta($fasta);
		
	# Loop over @matrix and get lenghts for the feature IDs in it
	MAT: for(my $cc=0; $cc<=$#matrix; $cc++)
		{
		my @arr = @{$matrix[$cc]};
		my $new_col = scalar(@arr);
		my $id = $matrix[$cc][$id_col];	

		# Loop over sequence lengths and find the length for the current geature ID
		LEN: for(my $dd=0; $dd<=$#length_matrix; $dd++)
			{
			if($id eq $length_matrix[$dd][0])
				{
				$matrix[$cc][$new_col] = $length_matrix[$dd][1];
				next MAT;
				}
			}
		} 

	# If the matrix had a header at the outset, add it back again
	if($header_y_n eq "y")	{	unshift(@matrix, \@header);	}
	
	if($file_or_var eq "file")	{	fileTools::write_table(\@matrix, "csv", "ln_${matrix_file_or_ref}", "lin");	}
	return(@matrix);
	}

return(1);

# end functions
