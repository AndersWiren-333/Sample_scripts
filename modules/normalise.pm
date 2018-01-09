package normalise;

# length_normalise_matrix($file_or_var, $filename_or_matrixref, $length_col, $header_y_n, $sample_col1, $sample_col2 ...)
# quantile_normalise_matrix($file_or_var, $filename_or_matrixref, $header_y_n, $tied_ranks_y_n, $preserve_zeros_y_n, $sample_col1, $sample_col2 ...)
# resampling_check()
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
require "$modules/normalise.pm";
require "$modules/listTools.pm";


# Create a timestamp string (can be attached to the name of logfiles, for example
my $timestamp = envir::timestamp();
my $rscript = "Rscript";

# end header

########################################################## Functions ##########################################################


sub length_normalise_matrix
	{
	# This function takes a gene expression matrix (either as a matrix reference – set $file_or_var to ‘var’ - or as a csv file – set $file_or_var to ‘file’)
	# and normalises expression values based on the length of genes. It returns the normalised matrix,  with expression values per 1000 bp of gene.
	# When running the function, specify the column number holding gene lengths ($length_col - if 0, set to "zero"), whether  the matrix has a header
	# ($header_y_n = "y", otherwise "n") and the column numbers holding the data to be normalised (i.e. those not containing annotation information).
	# The returned matrix contains all original headers and annotation information.

	# Set error messages and accept input parameters
	my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
	my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$file_or_var, \$filename_or_matrixref, \$length_col, \$header_y_n, \$sample_col1, \$sample_col2 ...)'\n\nwhere".
	"\t\$file_or_var should be set to 'var' if the input is a reference to a matrix, or 'file' if the input is a csv file (and the output also should be)\n".
	"\t\$filename_or_matrixref is either the name of the input csv file or a reference to the matrix to be normalised\n".
	"\t\$length_col is the number of the column holding length information (numbering starts at 1)\n".
	"\t\$header_y_n - set this to 'y' if the input file/matrix has a header or 'n' if it doesn't\n".
	"\t\$sample_col1, \$sample_col2 etc are the column numbers holding the data to be normalised (numbering starts at 1)\n\n";
	my @pars = @_ or die $usage;
	foreach my $el (@pars)  {       $el = text::trim($el);  }
	my $file_or_var = shift @pars or die $usage;
	my $matrix_file_or_ref = shift @pars or die $usage;      # NB! This will not work if the argument is the number 0 (because it will then be interpreted as false). In that case you need to use 'shift(@pars) or 0', but it may lead to proble$
	my $length_col = shift @pars or die $usage;
	my $header_y_n = shift @pars or die $usage;
	my @sample_cols = @pars;
	@sample_cols=misc::shift_input_cols(\@sample_cols);
	$length_col--;
	
	my @matrix=();
	if($file_or_var eq "file")	{	@matrix=fileTools::read_table($matrix_file_or_ref, "csv");	}
	elsif($file_or_var eq "var")	{	@matrix = @{$matrix_file_or_ref};	}
	else	{	die "\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Parameter \$file_or_var must be either 'file' or 'var'\n\n";	}
	if($length_col eq "zero")	{	$length_col = 0;	}
	
	# If there is a header, remove it from the matrix and store it in an array
	my @header=();
	if($header_y_n eq "y")
		{
		my $headerref = shift(@matrix);
		@header = @{$headerref};
		}

	# Loop over @matrix and normalise all values for each gene
	for(my $cc=0; $cc<=$#matrix; $cc++)
		{
		my $length = $matrix[$cc][$length_col];
		my $norm_factor = 1000/$length;		

		# Loop over sample expression values and normalise them
		for(my $dd=0; $dd<=$#sample_cols; $dd++)
			{
			my $col = $sample_cols[$dd];
			my $value = $matrix[$cc][$col];
			$matrix[$cc][$col] = $value * $norm_factor;
			}
		} 

	# If the matrix has a header at the outset, add it back again
	if($header_y_n eq "y")	{	unshift(@matrix, \@header);	}
	
	if($file_or_var eq "file")	{	fileTools::write_table(\@matrix, "csv", "lc_${matrix_file_or_ref}", "lin");	}
	return(@matrix);
	} # end length_normalise_matrix


sub quantile_normalise_matrix
	{
	# This function takes a gene expression matrix (either as a matrix reference – set $file_or_var to ‘var’ - or as a csv file –
	# set $file_or_var to ‘file’) and performs quantile normalisation on it (as described in Bolstad et al 2003 - but with additional
	# options). When running the function, specify whether the matrix has a header ($header_y_n = "y", otherwise "n") and the column
	# numbers holding the data to be normalised (i.e. those not containing annotation information. First column is called 1 (i.e.
	# not native perl indexing)). Make sure the column numbers are the last arguments to be specified. Quantile normalisation is a
	# rank based method, and you can choose whether to adjust tied ranks to the average rank within a tied group (set $tied_ranks to "y")
	# or not ($tied_ranks="n"). You can also choose whether zero-values in the original matrix should remain zero ($preserve_zeros = "y")
	# or be replaced by a normalised average expression value ($preserve_zeros = "n"). The returned matrix contains all original headers
	# and annotation information (if any). If the input matrix was read from a file, the normalised matrix will be written to an output
	# csv file (as well as returned as a matrix).
	
	# Set error messages and accept input parameters
	my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
	my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$file_or_var, \$filename_or_matrixref, \$header_y_n, \$tied_ranks_y_n, \$preserve_zeros_y_n, \$sample_col1, \$sample_col2 ...)'\n\nwhere".
	"\t\$file_or_var should be set to 'var' if the input is a reference to a matrix, or 'file' if the input is a csv file (and the output also should be)\n".
	"\t\$filename_or_matrixref is either the name of the input csv file or a reference to the matrix to be normalised\n".
	"\t\$header_y_n - set this to 'y' if the input file/matrix has a header or 'n' if it doesn't\n".
	"\t\$tied_ranks_y_n  If tied ranks should be replaced by average ranks within a tied group, set this to 'y', otherwise to 'n'.\n".
	"\t\$preserve_zeros_y_n  If values that are zero should not be replaced average values, set this to 'y', otherwise set it to 'n'\n".
	"\t\$sample_col1, \$sample_col2 etc are the column numbers holding the data to be normalised (first column is 1, second is 2 etc. i.e. not native perl indexing)\n\n";
	my @pars = @_ or die $usage;
	foreach my $el (@pars)  {       $el = text::trim($el);  }
	my $file_or_var = shift @pars or die $usage;
	my $matrix_file_or_ref = shift @pars or die $usage;      # NB! This will not work if the argument is the number 0 (because it will then be interpreted as false). In that case you need to use 'shift(@pars) or 0', but it may lead to proble$
	my $header_y_n = shift @pars or die $usage;
	my $tied_ranks = shift @pars or die $usage;
	my $preserve_zero = shift @pars or die $usage;
	my @sample_cols = @pars;
	my @shifted_sample_cols=misc::shift_input_cols(\@sample_cols);		# Translate column numbers to perl native indexing (where column 1 is column 0)

	print "\n\nReading matrix...  ";
	# Read/dereference the matrix
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
		}
	print "Finished!\n";
	
	print "Assigning ranks to and sorting column expression values...  ";
	# For each data column, get the ranks of genes from lowest to highest expression and record the ranks in @orig_rank_matrix,
	# sort the values in the column from lowest to highest and add the column to @sorted_expression_matrix
	my @orig_rank_matrix=();
	my @sorted_expression_matrix=();
	for(my $cc=0; $cc<=$#sample_cols; $cc++)
		{
		my $col = $sample_cols[$cc];
		my @data_column=matrixTools::get_matrix_columns_as_rows(\@matrix, $col);

		my @column_ranks=stats::assign_ranks(\@data_column, $tied_ranks, "num");
		@orig_rank_matrix=matrixTools::add_cols_to_matrix(\@orig_rank_matrix, \@column_ranks);
		
		@data_column = sort { $a <=> $b } @data_column;
		@sorted_expression_matrix=matrixTools::add_cols_to_matrix(\@sorted_expression_matrix, \@data_column);
		}
	print "Finished!\n";
	
	print "Computing means for rows in sorted expression matrix...   ";
	# Loop over rows in @sorted_expression_matrix and compute the (arithmetic) mean value of all columns in each row.
	# Store the row means in @row_means.
	my @row_means=matrixTools::apply(\@sorted_expression_matrix, "row", "mean");
	undef @sorted_expression_matrix;	# We no longer need this matrix, free some memory
	print "Finished!\n";
	
	print "Assigning ranks to row means...  ";
	# Get ranks for row means
	my @mean_ranks=stats::assign_ranks(\@row_means, $tied_ranks, "num");
	my @mat=();
	my @means_and_ranks_matrix=matrixTools::add_cols_to_matrix(\@mat, \@row_means, \@mean_ranks);
	@means_and_ranks_matrix = sort { $a->[1] <=> $b->[1] } @means_and_ranks_matrix;

	# Remove duplicates (if any) from @means_and_ranks_matrix (this is only to speed up processing - the full dataset will be returned in the end)
	my @u_means_and_ranks_matrix=misc::unique_matrix(\@means_and_ranks_matrix, 2, "num", "n");
	print "Finished!\n";
	
	# Go through @orig_rank_matrix and for every value (a rank for a gene expression value in a sample column), replace
	# the rank with the row mean that corresponds to that rank in @means_and_ranks_matrix.
	my @search_ranks = matrixTools::get_matrix_columns_as_rows(\@u_means_and_ranks_matrix, 2);
	
	print "Assigning normalised values to ranked original values...  \n\n";
	# Loop over rows in @orig_rank_matrix
	for(my $ee=0; $ee<=$#orig_rank_matrix; $ee++)
		{
		my $rad = 1+$ee;
		print "\tRow $rad ... ";
		my @row = @{$orig_rank_matrix[$ee]};
		
		# Loop over columns in row
		for(my $ff=0; $ff<=$#row; $ff++)
			{
			my $current_rank = $row[$ff];			
			my $search_rank_index=listTools::get_array_value_index(\@search_ranks, "num", $current_rank);
			
			my $new_value = "Valnot";
			
			# If the rank in @orig_rank_matrix isn't present in @search_ranks, get the two closest ranks
			if($search_rank_index eq "not_available")
				{
				my ($index1, $index2)=listTools::get_closest_indices(\@search_ranks, $current_rank);
				#print "Index1: $index1\tIndex2: $index2\n";
				my @low_high = ($u_means_and_ranks_matrix[${index1}-1][0], $u_means_and_ranks_matrix[${index2}-1][0]);
				#print "Taking means of: $low_high[0] , $low_high[1]\n";
				$new_value = stats::mean(\@low_high);
				}
			else	{	$new_value = $u_means_and_ranks_matrix[$search_rank_index][0];	}
			$orig_rank_matrix[$ee][$ff] = $new_value;
			}
		print "Done!\n";
		}
	print "\nFinished assigning normalised values to ranked original values\n";
	
	print "Setting those normalised values to 0 which should be...  ";
	# If zero-values in the original matrix should be preserved, i.e. not replaced by a row average (a normalised value greater than zero),
	# find the positions in the original matrix (@matrix) that have zeros. Set the corresponding positions in @orig_rank_matrix to zero.
	if($preserve_zero eq "y")
		{
		# Get relevant columns from original matrix
		my @selected_cols=matrixTools::get_matrix_columns(\@matrix, @sample_cols);

		my @dim_1=matrixTools::dim(\@selected_cols);
		my @dim_2=matrixTools::dim(\@orig_rank_matrix);
		my $equal_dims=listTools::array_equality(\@dim_1, \@dim_2);
		if($equal_dims eq "n")	{	die "The matrix holding new values and the matrix holding original values don't have the same dimensions\n\n\tSomething is wrong - try again\n\n";	}

		# Loop over rows in @selected_cols
		for(my $ff=0; $ff<=$#selected_cols; $ff++)
			{
			my @row2 = @{$selected_cols[$ff]};
			# Loop over columns in row
			for(my $gg=0; $gg<=$#row2; $gg++)
				{
				if($selected_cols[$ff][$gg] == 0)	{	$orig_rank_matrix[$ff][$gg] = 0;	}
				}
			}
		}
	print "Finished!\n";
	
	print "Checking that normalised matrix has the same dimensions as original matrix...  ";
	# Check that the normalised matrix has the proper dimensions and can be put back into the original matrix (which contains annotation information)
	my $num_cols_matrix=scalar(@shifted_sample_cols);
	my $num_rows_matrix=scalar(@matrix);
	my @dim_3 = ($num_rows_matrix, $num_cols_matrix);
	my @dim_4 = matrixTools::dim(\@orig_rank_matrix);
	my $equal_dims2=listTools::array_equality(\@dim_3, \@dim_4);
	if($equal_dims2 eq "n")	{	die "The matrix holding new values (corrected for zeros) and the matrix holding original values don't have the same dimensions\n\n\tSomething is wrong - try again\n\n";	}
	print "Finished!\n";
	
	print "Adding annotation information back to normalised matrix... ";
	# Replace the original columns in @matrix with @orig_rank_matrix
	my $firstcol = $sample_cols[0];
	my $lastcol = $sample_cols[$#shifted_sample_cols];
	my $dimstring = "1_${firstcol}_to_${num_rows_matrix}_${lastcol}";
	my @new_matrix=matrixTools::replace_matrix_chunk(\@matrix, $dimstring, \@orig_rank_matrix);
	print "Finished!\n";
	
	# If the matrix had a header at the outset, add it back again
	if($header_y_n eq "y")	{	unshift(@new_matrix, \@header);	print "Put header back onto matrix\n";	}
	
	print "Writing normalised matrix to outfile...  ";
	# If the matrix was read from an infile, print it to a new outfile
	if($file_or_var eq "file")	{	fileTools::write_table(\@new_matrix, "csv", "qn_${matrix_file_or_ref}", "lin");	}
	print "Finished!\n\n";
	print "\tAll done!\n\n";
	
	return(@new_matrix);
	} # end quantile_normalise_matrix

 
sub resampling_check
	{
	# Checks the effect on sample complexity and complexity of genome matching reads of resampling a set of non-redundant fasta files to a series of
	# consequtively smaller sizes (specified as percentages of original file size). Specify the size in percentage points between successive steps
	# (e.g. '5' means resample to 95%, 90%, 85% etc), the reference genome (fasta file) that reads should be aligned to, the number of mismatches to
	# allow in alignment between reference genome and an individual sequence read, and the names of the fasta files that should be checked.

	# Set error messages
	my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
	my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}).".
	" Correct usage: '${subname}(\$to_number, \$in_steps_of, \$align_to_this_genome_fasta, \$with_n_mismatches, \$fasta_file1, \$fasta_file2 etc...)\n\nwhere\n".
	"\t\$to_number is the number if reads to which files should be resampled\n".
	"\t\$in_steps_of is the numer of reads between each successive step of resampling of a particular fasta file\n".
	"\t\$align_to_this_genome_fasta is the genome fasta file to which resampled fasta files should be aligned to check number of genome matching reads.\n".
	"\t\$with_n_mismatches is the number of mismatches to allow in the alignment between genome and an individual read (if 0, say 'zero')\n".
	"\t\$fasta_file1, \$fasta_file2 etc... are the names of the files to be resampled\n\n";

	# Accept input parameters
	my @pars = @_ or die $usage;
	foreach my $el (@pars)  {       $el = text::trim($el);  }
	my $end_number = shift @pars or die $usage;
	my $steps_of = shift @pars or die $usage;
	my $ref_genome = shift @pars or die $usage;
	my $mismatch = shift @pars or die $usage;
	my @samples = @pars or die $usage;
	if($mismatch eq "zero")	{	$mismatch = 0;	}

	# Open a statistics outfile
	open(my $stat, ">>", "stats_subsampling_check.csv") or die "Subroutine '${subname}' (called by script '${calling_script}', ".
	"line ${calling_line}) couldn't create statistics outfile\n";
	my $stat_header = "File,Total_reads,Unique_reads,Complexity,Copies_per_read,Total_genome_matching_reads_${mismatch}mm,Perc_genome_matching_reads_${mismatch}mm,".
	"Unique_genome_matching_reads_${mismatch}mm,Perc_unique_genome_matching_reads,Complexity_genome_matching_reads_${mismatch}mm,Copies_per_read_genome_matching_reads_${mismatch}mm";
	print($stat "$stat_header\n");
	
	# Loop over samples
	for(my $c=0; $c<=$#samples; $c++)
		{
		my $infile = $samples[$c];
		my ($nr_fasta_gz, $nr_fasta)=misc::check_compressed($infile);
		
		# Set numbers to be resampled to
		my ($total_reads, $unique_reads)=fastaTools::count_nr_fasta_reads($nr_fasta);
		my $diff = ($total_reads-$end_number);
		my $nsteps = int($diff/$steps_of);	
		if($nsteps < 1)	{	next;	}

		my @numbers=();
		for(my $z=1; $z<=$nsteps; $z++)
			{
			my $increment = $z * $steps_of;
			my $number = $end_number+$increment;
			push(@numbers, $number);
			}
		push(@numbers, $end_number);
		@numbers = sort { $b <=> $a } @numbers;
	

		# Subsample to a range of numbers (loop over @numbers)
		for(my $d=0; $d<=$#numbers; $d++)
			{
			my @stats_vector=();			# Set up an array to hold statistics for this particular subsampled file
			my $numA = $numbers[$d];
			print "$nr_fasta\t$numA\n";
			my $r_sub_fasta = "r_${numA}_${nr_fasta}";

			my $patfile = "frigg.pat";
			#$patfile =~ 's/fa/pat/';

			# Subsample
			system("perl $scripts/bootstrap.pl $nr_fasta -1 $numA 1 $r_sub_fasta");
			system("perl $scripts/R2NR.pl $r_sub_fasta frigg.fa");

			# Count reads in resulting file
			my ($total_reads, $uni_reads)=fastaTools::count_nr_fasta_reads("frigg.fa");
			my $complexity=($uni_reads/$total_reads);
			my $cop_per_read=1/$complexity;
			push(@stats_vector, $r_sub_fasta, $total_reads, $uni_reads, $complexity, $cop_per_read);

			# Align to genome
			system("patman -D $ref_genome -P frigg.fa -o $patfile -e $mismatch");
			unlink($r_sub_fasta);
			unlink("frigg.fa");

			# Count genome matching reads
			my ($gen_total_reads, $gen_uni_reads, $gen_complexity, $gen_cop_per_read, $countM)=misc::count_aligned_reads($patfile, "nr");

			my $perc_gen_total_reads = ($gen_total_reads/$total_reads)*100;
			my $perc_gen_uni_reads = ($gen_uni_reads/$uni_reads)*100;
			push(@stats_vector, $gen_total_reads, $perc_gen_total_reads, $gen_uni_reads, $perc_gen_uni_reads, $gen_complexity, $gen_cop_per_read);
			my $out_stats=join(",", @stats_vector);
			print($stat "$out_stats\n");
			unlink($patfile);

			} # Loop over percentages ends here

		system("gzip $nr_fasta");
		} # Loop over samples ends here

	close($stat);
	} # end resampling_check


return(1);

# end functions
