package matrixTools;

# sub testprint_matrix(\@matrix)
# sub transpose_matrix(\@matrix)
# sub generate_rand_matrix($rows, $columns, $min, $max, $decimals)
# sub create_zero_matrix($num_rows, $num_cols)
# sub merge_matrices_col(\@matrix1)
# sub split_matrix_on_col_value(\@matrix, $col_index)
# sub get_matrix_column($matrixref, $col_index1, $col_index2...)
# sub delete_matrix_columns($matrixref, $col_index1, $col_index2...)
# sub add_length_to_matrix_gene_id($matrixref, $length_col_index, $header_y_n)
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

# Set paths to scripts and modules. Setting explicit paths to the scripts and modules in this specific repository (rather than adding paths to @INC, PERLLIB and PATH on
# your local system) avoids the risk of scripts calling the wrong scripts/modules if you have other repositories on your system that happen to have some script- and module names
# in common with this repository.
my $scripts = cwd;
chdir("../modules");
my $modules = cwd;
chdir("../scripts");

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

# Prints a given matrix on screen
sub testprint_matrix
        {
	my $usage = "Syntax error for sub testprint_matrix. Correct usage: 'matrixTools::testprint_matrix(\@matrix)'\n";
        my $matrix_ref = $_[0] or die $usage;
        my @matrix = @{$matrix_ref};

        # Loop over rows in matrix
        for(my $c=0; $c<=$#matrix; $c++)
                {
                my @arr=@{$matrix[$c]};
		my $last_col = scalar(@arr);		
		for(my $d=0; $d<$last_col; $d++)
			{
			print("$matrix[$c][$d]\t");
			}

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
# a specific column. NB! The first column in a matrix has to be specified as "n" rather than "0", since Perl can't accept "0" as a command line argument.
sub get_matrix_columns
        {
        # Set error messages and accept input parameters
        my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
        my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$matrixref, \$col_index1, \$col_index2 etc)'\n\nwhere".
        "\t\$matrixref is a refernce to the matrix from which columns should be retrieved\n".
        "\t\$col_index1 is the index of the first column to be retrieved\n".
	"\t\$col_index2 is the index of the second column to be retrieved, etc.\n".
	"\n".
	"\tNB! The first column in a matrix has to be specified as 'n' rather than '0', since Perl can't accept '0' as a command line argument.\n\n";
        my @pars = @_ or die $usage;
        foreach my $el (@pars)  {       $el = text::trim($el);  }
        my $matrixref = shift @pars or die $usage;      # NB! This will not work if the argument is the number 0 (because it will then be interpreted as false). In that case you need to use 'shift(@pars) or 0', but it may lead to proble$
        my @cols = @pars or die $usage;
	my $num_cols = scalar(@cols);
	my @matrix = @{$matrixref};

	# If @cols contains "n" (the represenation of zero), replace it with 0.
	foreach my $el (@cols)
		{
		if($el eq "n")	{	$el = 0;	}
		}

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

# Deletes specified columns from a matrix. NB! The first column in a matrix has to be specified as "n" rather than "0", since Perl can't accept "0" as a command line argument.
sub delete_matrix_columns
        {
        # Set error messages and accept input parameters
        my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
        my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$matrixref, \$col_index1, \$col_index2 etc)'\n\nwhere".
        "\t\$matrixref is a reference to the matrix from which columns should be deleted\n".
        "\t\$col_index1 is the index of the first column to be deleted\n".
        "\t\$col_index2 is the index of the second column to be deleted, etc.\n".
        "\n".
        "\tNB! The first column in a matrix has to be specified as 'n' rather than '0', since Perl can't accept '0' as a command line argument.\n\n";

        my @pars = @_ or die $usage;
        foreach my $el (@pars)  {       $el = text::trim($el);  }
        my $matrixref = shift @pars or die $usage;      # NB! This will not work if the argument is the number 0 (because it will then be interpreted as false). In that case you need to use 'shift(@pars) or 0', but it may lead to proble$
        my @exclude_col_indices = @pars or die $usage;
        my @matrix = @{$matrixref};
	
        # If @exclude_col_indices contains "n" (the represenation of zero), replace it with 0.
        foreach my $el (@exclude_col_indices)
                {
                if($el eq "n")  {       $el=0;       }
                }

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

return(1);

# end functions
