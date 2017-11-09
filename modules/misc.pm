package misc;

# sub patman_align($fasta_infile, $patman_outfile_name, $ref_genome, $max_seqs_per_file, $mismatches)
# sub read_filelist($filelist)
# sub file_to_matrix($file, $separator, $num_rows, $col_1_name, $col_2_name etc)
# sub unique_matrix(\@matrix)
# sub unique_list(\@array)
# sub count_aligned_reads($patman_output_file)
# sub blast_to_matrix($blastfile.tsv, $top_num_hits, $exclude_string1, $exclude_string2 etc...)
# sub check_compressed($filename)
# sub length_normalise_matrix(\@matrix)
# sub length_from_id_to_column($infile, $outname)
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

# Create a timestamp string (can be attached to the name of logfiles, for example
my $timestamp = envir::timestamp();
my $rscript = "Rscript";

# end header

########################################################## Functions ##########################################################

# Aligns DNA or RNA sequences in a fasta file to a reference sequence (e.g. genome or transcriptome) allowing a given number of mismatches between sequence and reference. 
# Automaticllay splits infiles into smaller pieces if they are larger than a given number of sequences. The output is an alignment information file (in patman format). The
# function then counts the total and unique (non-redundant) number of reference matching reads, the complexity and number of copies per read and the number of reads
# that map to the reference sequence in multiple places, and returns these numbers as an array.
sub patman_align
	{
        # Set error messages and accept input parameters
        my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
        my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: ".
	"'${subname}(\$fasta_infile, \$patman_outfile_name, \$ref_genome, \$max_seqs_per_file, \$mismatches, \$fasta_format)'\n\nwhere".
        "\t\$fasta_infile is the fasta file with reads to be aligned\n".
        "\t\$patman_outfile_name is the desired name of the outfile\n".
        "\t\$ref_genome is the reference sequence to align to\n".
	"\t\$max_seqs_per_file - fasta infiles with more sequences will be split into pieces before alignment to conserve memory\n".
	"\t\$mismatches is the maximum allowed differences between a sequence read and the reference sequence\n".
	"\t\$fasta_format is the format of the fasta files to be aligned. This can be either 'r', 'nr' or 'mnr'\n\nwhere".
	"\t'r' is 'redundant' i.e. '>id (linebreak) sequence'\n".
	"\t'nr' is 'non-redundant' i.e. '>sequence-abundance (linebreak) sequence'\n".
	"\t'mnr' is 'multi-non-redundant' i.e. '>sequence_abundance1_abundance2_abundance3... (linebreak) sequence'\n\n";
        my @pars = @_ or die $usage;
        #foreach my $el (@pars)  {       $el = text::trim($el);  }

	my $infile_fa = shift @pars or die $usage;
	my $outfile_pat = shift @pars or die $usage;
	my $ref_genome = shift @pars or die $usage;
	my $seqs_per_piece = shift @pars or die $usage;
	my $mismatch = shift @pars or 0;
	my $fasta_format = shift(@pars) or die $usage;
	my $temp_pat = "$outfile_pat"."_TEMP";

        my $num_parts = fastaTools::divide_fasta($infile_fa, $seqs_per_piece);

        for(my $d=1; $d<=$num_parts; $d++)
                {
                my $partfile_fa = "${infile_fa}_${d}.fa";
                my $partfile_pat = "$partfile_fa".".pat";
                system("patman -D $ref_genome -P $partfile_fa -o $partfile_pat -e $mismatch");
                system("cat $partfile_pat >> $temp_pat");
                unlink($partfile_fa, $partfile_pat);
                }
	move($temp_pat, $outfile_pat);
	my @stats=misc::count_aligned_reads($outfile_pat, $fasta_format);	# This assumes that the input fasta files were in non-redundant format, doesn't it?
	return(@stats);
	}

# Reads a given plain text file with filenames (one on each line) into a matrix.
sub read_list
	{
	# Set error messages and accept input parameters
        my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
	my $usage="\nSyntax error for sub ${subname}. Correct usage: '${subname}(\$filelist)'\n";
	my $infile = $_[0] or die $usage;
	if(!open(IN, $infile))	{	die "sub $subname couldn't open filelist $infile\n";	}
	my @names=();

	while(<IN>)
		{
		my $file=text::trim($_);
		push(@names, $file);
		}	
	return(@names);
	}

# Reads the given number of rows and the named columns of a given file (with a given type of field separator) into a matrix. The return value is a matrix with three lines.
# The first line is a reference to the whole matrix (without headers), the second line is a reference to an array holding the headers (if the original file had no headers,
# this array will hold headers of the form "column_0, column_1...") of the selected columns and the third line is a reference to an array holding the indices (column numbers)
# of the selected columns.
sub file_to_matrix
        {
        # Set error messages and accept input parameters
        my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
        my $usage = "\nSyntax error for sub $subname. Correct usage: '${subname}(\$file, \$separator, \$headers_y_n, \numrows=333, \$col_1_name, \$col_2_name etc)'\n\nwhere\n".
	"\t\$file is the infile (required)\n".
	"\t\$separator is which type of delimiter used between fields in the infile (options are 'tab', 'comma' and 'semicolon')(required)\n".
	"\t\$headers_y_n is an indicator of whether the input matrix has headers or not (options are 'y' for yes and and 'n' for no)(required)\n".
	"\t\$num_rows is the number of rows to read from the chosen file. (optional, default=all rows)\n".
	"\t\$col_1_name is the name of the first column you want to read (if the infile has headers), or the column number of the first column to read (if the infile doesn't have headers)(optional, default=all columns)\n".
	"\t\$col_1_name is the name or number of the first column you want to read (optional, default=all columns)\n\n".
	"\tThe subroutine returns a matrix where the first line is a reference to the whole matrix (without headers), the second line is a reference to an array holding the headers\n".
	"\t(if the infile has no headers, the returned headers will be 'column_1, column_2' etc.) and the third line is a reference to an array holding the indices (column numbers) of the desired columns.\n\n";

	# Accept required input parameters
        my @pars = @_;
	foreach my $el (@pars)	{	$el = text::trim($el);	}
        my $infile = shift(@pars) or die $usage;
	my $sep = shift(@pars) or die $usage;
	my $headers_y_n = shift(@pars) or die $usage;
	unless(($sep eq "tab") or ($sep eq "comma") or ($sep eq "semicolon"))  {       die "Syntax error for sub $subname. Argument \$separator must be either 'tab', 'comma' or 'semicolon'. Try again!\n";    }
	unless(($headers_y_n eq "y") or ($headers_y_n eq "n"))  {       die "Syntax error for sub $subname. Argument \$headers_y_n must be either 'y' or 'n'. Try again!\n";    }
	my $nrows="all";
	my $all_cols_y_n = "y";

	# Set $nrows if specified (else it remains "all")
	if((defined($pars[0])) and ($pars[0] =~ m/numrows/))
		{
		$nrows=shift(@pars);
		my @rowparts = split("=", $nrows);
		$nrows = $rowparts[1];
		}

        my @sel_colnames=();
        my @sel_colindices=();

	# Set column names/indices if specified (else $all_cols_y_n remains "y")
	if(defined($pars[0]))
		{
		if($headers_y_n eq "y")	{	@sel_colnames = @pars;	}
		if($headers_y_n eq "n") {       @sel_colindices = @pars;  }		
		$all_cols_y_n = "n";
		}


	my @headers=();
        my $linenum=0;

	if($sep eq "tab")	{	$sep = "\t";	}
	elsif($sep eq "comma")	{	$sep = ",";	}
	elsif($sep eq "semicolon")  {       $sep = ";";     }

	# Open infile
        if(!open(IN, $infile))  {       die "sub $subname couldn't open infile $infile\n"; }
        my @matrix_all_cols=();
        my $lineno=0;
	my $stop="alfa";
	if($nrows ne "all")
		{	
		$stop = $nrows+1;
		if($headers_y_n eq "y") {       $stop++;        }
		}

        # Loop over lines in file
	my @colnames=();
        LINES: while(<IN>)
                {
		$lineno++;
                my $line=text::trim($_);

		# Make sure that empty lines are not included in the matrix (these may contain formatting characters and therefore
		# be interpreted as data though they are empty if you are transferring data from an excel sheet)
		my $empty="y";
		if($line =~ m/\w/)  {       $empty="n";     }
		if($empty eq "y")	{	$lineno--; next LINES;	}

		# Make sure to skip comment lines
		if($line =~ m/^#/m)	{	$lineno--; next LINES;	}

		# Process current line
                my @arr = split("$sep", $line);
                foreach my $le (@arr)   {       $le=text::trim($le);    }

		# If we have passed the number of rows that should be returned, exit loop at this stage
                if(($nrows ne "all") and ($lineno==$stop))   {       last;   }

                # If there is a header row, read it
                if(($headers_y_n eq "y") and ($lineno==1))	{	@colnames = @arr; next;	}

		# If not, make a header row
		elsif(($headers_y_n eq "n") and ($lineno==1))
			{
			for(my $cc=0; $cc<=$#arr; $cc++)
				{
				my $num = $cc;
				my $name = "column_".$num;
				push(@colnames, $name);
				}
			}
		push(@matrix_all_cols, \@arr);
                } # Loop over lines in file ends here
	close(IN);

	# Now we have the whole matrix and an array with column names. 
	# If only selected columns are to be returned
		# and there are headers, we have the names of the columns but not the indices
		# and there are no headers, we have the indices of the columns, and constructed names of them

	my @matrix=();

	# If all columns should be returned
	if($all_cols_y_n eq "y")
		{
		@sel_colnames = @colnames;
		@sel_colindices = (0..$#sel_colnames);
		@matrix = @matrix_all_cols;
		@matrix_all_cols=();
		}

	# Else, if only selected columns should be returned
	elsif($all_cols_y_n eq "n")
		{
		# Get indices of the selected column names
		if($headers_y_n eq "y")
			{
			for(my $dd=0; $dd<=$#sel_colnames; $dd++)
				{
				my $name = $sel_colnames[$dd];
				($sel_colindices[$dd]) = grep { $colnames[$_] =~ m/$name/ } 0..$#colnames;
				}
			}

		# Or get names of the selected column indices	
		elsif($headers_y_n eq "n")
			{
			for(my $ee=0; $ee<=$#sel_colindices; $ee++)
                                {
                                my $num = $sel_colindices[$ee];
                                $sel_colnames[$ee] = "column_".$num;
                                }
			}

		
		# Loop over @matrix_all_cols and send the selected columns to @matrix
		for(my $ff=0; $ff<=$#matrix_all_cols; $ff++)
			{
			my @arr = @{$matrix_all_cols[$ff]};
			my @new_arr=();
			for(my $gg=0; $gg<=$#sel_colindices; $gg++)
				{
				my $index = $sel_colindices[$gg];
				push(@new_arr, $arr[$index]);
				}
			push(@matrix, \@new_arr);
			}
		@matrix_all_cols=();
		} # End of "if only selected columns should be returned"-statement

	my @grand_matrix = (\@matrix, \@sel_colnames, \@sel_colindices);
        return(@grand_matrix);
        }


# Removes duplicates from a matrix
sub unique_matrix
	{
        # Set error messages and accept input parameters
        my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
	my $usage = "\nSyntax error for sub ${subname}. Correct usage: '${subname}(\$matrixref, \$column, \$num_or_alph, \$headers_y_n)'\n\nwhere\n".
        "\t\$matrixref is a reference to the matrix you want to check for duplicates\n".
	"\t\$column is the index of the column you want to base the uniqueness on (where column 1 has index 1, 2 has index 2 etc.\n)".
	"\t\$num_or_alph indicates whether the values in the chosen column are numerical or alphabetic ('num' or 'alph')\n".
	"\t\$headers_y_n is an indicator of whether the input matrix has headers or not. Options are 'y' for yes and 'n' for no.\n\n";

	# Accept input and parameters
	my @pars = @_;
	my $matrix_ref = text::trim($pars[0]);
	my $col = text::trim($pars[1]);
	$col--;	# Make sure column index is in perl style (in perl, the first column is always called 0, the second 1 etc.)
	my $num_or_alph = text::trim($pars[2]);
	my $headers_y_n = text::trim($pars[3]);
	if(($num_or_alph ne "num") and ($num_or_alph ne "alph"))	{	die "The \$num_or_alph argument for sub $subname must be either 'num' or 'alph'. Try again!\n";	}
	my @headers=();

	# Read and sort the matrix
	my @matrix = @{$matrix_ref};	
	if($headers_y_n eq "y")
		{
		my $header_ref = shift(@matrix);
		@headers = @{$header_ref};
		}
	if($num_or_alph =~ m/num/)	{	@matrix = sort { $a->[$col] <=> $b->[$col] } @matrix;	}
	elsif($num_or_alph =~ m/alph/)	{	@matrix = sort { $a->[$col] cmp $b->[$col] } @matrix;	}

	# Loop over lines in matrix
	my @new_matrix=();
	my $comp_val = $matrix[0][$col];
	push(@new_matrix, $matrix[0]);

	for(my $cc=1; $cc<=$#matrix; $cc++)
		{
		my $value = $matrix[$cc][$col];
		if($num_or_alph =~ m/num/)
			{
			# If new value, add to @new_matrix
			unless($value==$comp_val)
				{
				push(@new_matrix, $matrix[$cc]);
				$comp_val=$value;
				}
			}
		if($num_or_alph =~ m/alph/)
                        {
			# If new value, add to @new_matrix
                        unless($value eq $comp_val)
                                {
                                push(@new_matrix, $matrix[$cc]);
                                $comp_val=$value;
                                }
                        }
		}
	return(@new_matrix);
	}


# Removes duplicates from an array
sub unique_list
        {
        # Set error messages and accept input parameters
        my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
        my $usage = "\nSyntax error for sub ${subname}. Correct usage: '${subname}(\$arrayref, \$num_or_alph)'\n\nwhere\n".
        "\t\$arrayref is a reference to an aray or matrix\n".
        "\t\$num_or_alph indicates whether the values in the chosen column are numerical or alphabetic ('num' or 'alph')\n\n";
        my @pars = @_ or die $usage;
	foreach my $el (@pars)	{	$el = text::trim($el);	}
        my $arref = $pars[0];
        my $num_or_alph = $pars[1];
        if(($num_or_alph ne "num") and ($num_or_alph ne "alph"))        {       die "The \$num_or_alph argument for sub $subname must be either 'num' or 'alph'. Try again!\n";      }

        my @arr = @{$arref};
        if($num_or_alph =~ m/num/)      {       @arr = sort { $a <=> $b } @arr;   }
        elsif($num_or_alph =~ m/alph/)  {       @arr = sort { $a cmp $b } @arr;   }

        # Loop over lines in matrix/arr
        my @new_arr=();
        my $comp_val = $arr[0];
        push(@new_arr, $arr[0]);

        for(my $cc=1; $cc<=$#arr; $cc++)
                {
                my $value = $arr[$cc];
                if($num_or_alph =~ m/num/)
                        {
                        unless($value==$comp_val)
                                {
                                push(@new_arr, $arr[$cc]);
                                $comp_val=$value;
                                }
                        }
                if($num_or_alph =~ m/alph/)
                        {
                        unless($value eq $comp_val)
                                {
                                push(@new_arr, $arr[$cc]);
                                $comp_val=$value;
                                }
                        }
                }
        return(@new_arr);
        }

# Counts the number of total, unique and multiply aligning reads in a patman alignment file and returns those numbers along with sample complexity and
# copies per read. This is the array returned: @($countT, $countU, $complex, $cop_per_read, $countM)
sub count_aligned_reads
	{
        # Set error messages and accept input parameters
        my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
	my $usage = "\nUsage error for sub ${subname} (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$patman_output_file, \$fasta_format)'\n\nwhere\n".
        "\t\$patman_output_file is a file containing information about where sequence reads align to a reference sequence\n".
	"\t\$fasta_format is the format of the fasta files to be aligned. This can be either 'r', 'nr' or 'mnr'\n\nwhere".
	"\t\t'r' is 'redundant' i.e. '>id (linebreak) sequence'\n".
	"\t\t'nr' is 'non-redundant' i.e. '>sequence-abundance (linebreak) sequence'\n".
	"\t\t'mnr' is 'multi-non-redundant' i.e. '>sequence_abundance1_abundance2_abundance3... (linebreak) sequence'\n\n";
        my @pars = @_ or die $usage;
        foreach my $el (@pars)  {       $el = text::trim($el);  }
        my $infile = shift @pars or die $usage;
	my $fasta_format = shift @pars or die $usage;

	open(my $in, "<", $infile) or die "Subroutine $subname (called by script '${calling_script}', line ${calling_line}) couldn't open infile $infile\n";
	my @matrix=();
	while(my $line = <$in>)
		{
		$line = text::trim($line);
		my @arr = split(/\t/, $line);
		#foreach my $el (@arr)	{	$el=text::trim($el);	}	# I get a usage error for text::trim here fore no apparent reason (does text::trim require $el to be a string rather than a number?)
		push(@matrix, \@arr);		
		}
	close($in);

#	my ($matrix_ref, $header_ref, $col_indices_ref)=misc::file_to_matrix($infile, "tab", "n");
#	my @matrix=@{$matrix_ref};

	# Sort the patman file matrix on ID
	@matrix = sort { $a->[1] cmp $b->[1] } @matrix;

        my $countT=0;
        my $countU=0;
        my $countM=0;
	my $multiple="n";
        my $complex = 0;
        my $cop_per_read = 0;

	if(@matrix)
		{
		if($fasta_format eq "r")
			{
			my @comparr=@{$matrix[0]};
			if(defined($comparr[1]))
				{
				$countT++;
                		$countU="unknown";
			
				# Loop over reads
				R_COUNT: for(my $cc=1; $cc<=$#matrix; $cc++)
					{
					if($matrix[$cc][1] eq $comparr[1])
						{
						$multiple="y";
						if($cc == $#matrix)	{	$countM++;	}
						next R_COUNT;				
						}
					else
						{
						$countT++;
						if($multiple eq "y")    {       $countM++;      }
						$multiple="n";
						@comparr=@{$matrix[$cc]};
						}
					} # Loop over reads ends here
				} 
			} # If r format ends here

        	elsif($fasta_format eq "nr")
                	{
			my @comparr=@{$matrix[0]};
			if(defined($comparr[1]))
				{
				$countU++;
				if($comparr[1] =~ m/-(\d+)/)    {       $countT = $1;   }

				# Loop over reads
				NR_COUNT: for(my $cc=1; $cc<=$#matrix; $cc++)
					{
	                        	if($matrix[$cc][1] eq $comparr[1])
        	                        	{
                	                	$multiple="y";
						if($cc == $#matrix)     {       $countM++;      }
                        	        	next NR_COUNT;
                                		}
                        		else
                                		{
                                		$countU++;
                                		if($multiple eq "y")    {       $countM++;      }
                                		$multiple="n";
                                		if($matrix[$cc][1] =~ m/-(\d+)/)
                                        		{
                                        		my $num = $1;
                                        		$countT = $countT+$num;
                                        		}
						@comparr=@{$matrix[$cc]};
                                		}
					} # Loop over reads ends here
				}
			} # If nr format ends here

        	elsif($fasta_format eq "mnr")
                	{
                	my @comparr=@{$matrix[0]};
                	if(defined($comparr[1]))
                        	{
                        	$countU++;
				my @comp_abundances=split("_", $comparr[1]);
				shift(@comp_abundances);
				foreach my $el (@comp_abundances)
					{
					#$el = text::trim($el);
					$countT = $countT+$el;
					}

                        	# Loop over reads
                        	MNR_COUNT: for(my $cc=1; $cc<=$#matrix; $cc++)
                                	{
                                	if($matrix[$cc][1] eq $comparr[1])
                                        	{
                                        	$multiple="y";
                          	              	if($cc == $#matrix)     {       $countM++;      }
                                	        next MNR_COUNT;
                                        	}
                    	            	else
                        	                {
                                	        $countU++;
                                        	if($multiple eq "y")    {       $countM++;      }
                                        	$multiple="n";

                 		       		my @abundances=split("_", $matrix[$cc][1]);
						shift(@abundances);
                        			foreach my $el (@abundances)
                                			{
                                			#$el = text::trim($el);		# Usage error for text::trim again here		
                                			$countT = $countT+$el;
                                			}						
						@comparr=@{$matrix[$cc]};
                                        	}
                                	} # Loop over reads ends here
                        	}
			} # If mnr format ends here

		if($countU eq "unknown")
			{
			$complex = "unknown";
			$cop_per_read = "unknown";
			}
		else
			{
			unless($countT==0)	{	$complex = $countU/$countT;	}
			unless($countU==0)	{	$cop_per_read = $countT/$countU;	}
			}
	
		} # End of if defined @matrix
	
	else	{	($countT, $countU, $complex, $cop_per_read, $countM) = (0,0,0,0,0);	}

	my @results = ($countT, $countU, $complex, $cop_per_read, $countM);
	return(@results);
	}

# Reads a textfile containing blasthits into a matrix (alphabetically sorted on sequence id, ascending) and returns the matrix. One row for each blasted sequence, first column is sequence id, all other
# columns are blast hits in order of highest bitscore, % identity, alignment length lowest e-value 
sub blast_to_matrix
	{
        # Set error messages and accept input parameters
        my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
        my $usage="\nSyntax error for sub ${subname}. Correct usage: '${subname}(\$blastfile.tsv, \$top_num_hits, \$exclude_string1, \$exclude_string2 etc...)'\n\nwhere".
        "\t\$blastfile.tsv is a file with blast results, in tsv format\n".
        "\t\$top_num_hits is the number of best hits to report for each blasted sequence\n".
        "\t\$exclude_string1 (etc) are any number of (optional) strings that, if a blast hit contains them in its name, excludes that blast hit from being reported\n\n";
        my @pars = @_ or die $usage;
        foreach my $el (@pars)  {       $el = text::trim($el);  }
        my $blastfile = shift @pars or die $usage;
	my $topnum = shift @pars or die $usage;
	my @exclude_strings=();
	if(defined($pars[0]))	{	@exclude_strings=@pars;	}

	# Open filehandle(s)
	open(my $blast, "<", $blastfile) or die "sub $subname couldn't open infile $blastfile\n";

	# Read blast result file
	my @blast_hits=();
	while(my $line = <$blast>)
        	{
        	$line = text::trim($line);
        	my @arr=split("\t", $line);
        	for(my $i=0; $i<=$#arr; $i++)    {       $arr[$i]=text::trim($arr[$i]);        }
        	$arr[1] =~ s/,//g;	# If the blast hit name contains comma, remove it, because it can cause a formatting error if the results are to be written to a csv file later (and they probably are) 
        	$arr[0] =~ s/u//;	# If the Feature_ID of the blast hit has an initial "u" in it (in some parts of DE ananlysis it is used to indicate that this is an unannotated transcript), remove the "u". 
        	push(@blast_hits, \@arr);
        	}
	close($blast);

	# Sort blast hits on Feature_ID, bitscore, % identity, alignment length and e-value
	@blast_hits = sort { $a->[0] cmp $b->[0] || $b->[3] <=> $a->[3] || $b->[6] <=> $a->[6] || $b->[5] <=> $a->[5] || $a->[2] <=> $b->[2] } @blast_hits;

	# The format of a line in @blast_hits is now:

	#       0       Query seqid             (The id line of your blasted fasta sequence)					(first sorting criterion)
	#       1       Subject title		(The description of the blast hit)
	#       2       E-value                 (Indication of how likely it is to find a hit for that sequence by chance)	(fifth sorting criterion)
	#       3       Bitscore                (Adjusted score of how good the hit is)						(second sorting criterion)
	#       4       Raw score		(Raw score of how good the hit is)
	#       5       Alignment length	(Number of basepairs where query and hit align)        				(fourth sorting criterion)
	#       6       % identity              (Percentage of aligned basepairs where query and hit are identical)		(third sorting criterion)


	# If textstrings have been specified that should exclude a blast hit from being reported if it contains them, remove any hits
	# with those strings from the matrix.
	my @new_blast_hits=();
	if(defined($exclude_strings[0]))
		{
		# Loop over lines in @blast_hits
		for(my $cc=0; $cc<=$#blast_hits; $cc++)
        		{
			my $keep="y";
			
			# Loop over @exclude_strings
			for(my $dd=0; $dd<=$#exclude_strings; $dd++)
				{
				if($blast_hits[$cc][1] =~ m/$exclude_strings[$dd]/)	{	$keep="n";	}
				}

			# Add to @new_blast_hits
			if($keep eq "y")	{	push(@new_blast_hits, $blast_hits[$cc]);	}
        		}
		@blast_hits=@new_blast_hits;
		@new_blast_hits=();
		}

	# Reformat the blast hit matrix so that each gene has one entry, and the format of that entry is this:

	#       0       	Gene_id
	#       1       	Name of blast-hit 1
	#       2       	Name of blast-hit 2
	#	..		..
	#       $topnum		Name of blast-hit $top_num


	my @new_blast_matrix=();        # Make new matrix to hold the formatted list
	my $comp_id="";
	my @outarr=();  # Initialize a small array to hold a gene id and that genes' top $topnum blast hits

	# Loop over @blast_hits
	for(my $ee=0; $ee<=$#blast_hits; $ee++)
        	{
	        my @arr=@{$blast_hits[$ee]};
        	my $id=$arr[0];

		# If this is the first line of the hit matrix...
	        if($ee==0)
        	        {
                	$comp_id=$id;
                	@outarr=($arr[0], $arr[1]);             # Start outline with Feature_ID and hit name (since this is the first (=top) hit for that gene)
                	}

		# If this is any other line...
        	else
                	{
			# If we are still on the same gene as previously...
                	if($id eq $comp_id)
                        	{
                        	push(@outarr, $arr[1]);         # Add the name of the current hit to the growing hit list for the current gene
                        	if($ee==$#blast_hits)                    # If this is the last line of the hit matrix...
                                	{
					my $new_topnum=$topnum;
					if($#outarr < $topnum)	{	$new_topnum = $#outarr;	}
                                	my @newarr=@outarr[0..$new_topnum];               # ...add the top $topnum blast hits for the current gene to the new, reformatted blast matrix
                                	push(@new_blast_matrix, \@newarr);
                                	}
                        	}

			# If we are on a new gene...
                	else
                        	{
				my $new_topnum=$topnum;
				if($#outarr < $topnum)  {       $new_topnum = $#outarr;	}
                        	my @newarr=@outarr[0..$new_topnum];               # add the top $topnum blast hits for the previous gene to the new, reformatted blast matrix
                        	push(@new_blast_matrix, \@newarr);
                        	$comp_id=$id;                   # ...start growing the hit list for the new gene.
                        	@outarr=($arr[0], $arr[1]);
                        	if($ee==$#blast_hits)
                                	{
                                	push(@new_blast_matrix, \@outarr);	# We use @outarr here because if we are on the last line and this is a new gene, then it it follows that the gene only has one blast hit.
                                	}
                        	}

                	}
        	} # Loop over @blast_hits ends here
	@blast_hits=();

	# Print matrix to csv file
	matrixTools::matrix_to_csv(\@new_blast_matrix, "top_${topnum}_blast_hits_per_gene.csv");

	# Return the matrix
	return(@new_blast_matrix);	
	}

# Checks whether a file is compressed and if so uncompresses it and returns the name of the compressed file
# as well as the uncompressed file. 
sub check_compressed
	{
        # Set error messages and accept input parameters
        my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
        my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$infile)'\n\nwhere".
        "\t\$infile is the file to be checked\n\n";
        my @pars = @_ or die $usage;
        foreach my $el (@pars)  {       $el = text::trim($el);  }
        my $infile = shift @pars or die $usage;

	my ($dir, $file, $basename, $suffix)=text::parse_path($infile);

	my $filename_gz="";
	my $filename="";

	# If the filename referes to a compressed file...
	if($suffix eq "gz")
		{
		# If a file with this name exists...
		if(-e $infile)
			{
			$filename_gz = $file;	# Set compressed name
			$filename = $basename;	# Set uncompressed name
			gunzip "$filename_gz" => "$filename"; unlink($filename_gz);		# Uncompress
			}

		# If it doesn't, check if there is an uncompressed file with the same basename
		else
			{
			# If there such an uncompressed file
			if(-e $basename)
				{
				$filename_gz = $file;
				$filename = $basename;
				}
			# If there isn't
			else	{	die "Subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}) couldn't find infile $infile\n";	}
			}
		}

	# If the filename referes to an uncompressed file...
	else
		{
		$filename = $file;
		$filename_gz = "$filename".".gz";

		# If there is no file with this name, check if there is a compressed file with the same basename
		unless(-e $infile)
			{
			# If there is such a file, uncompress it
			if(-e $filename_gz)	{	gunzip "$filename_gz" => "$filename"; unlink($filename_gz);	}
			}

		}

	my @outarr=($filename_gz, $filename);
        return(@outarr);
	}
	
# This function takes a gene expression matrix and normalises expression values based on the length of genes. It returns the normalised matrix, with expression values per 1000 bp of gene.
# The matrix must have a header, contain a column called "length", the first column must be the gene identifier, and the matrix must not contain any other columns than these two and the
# expression value columns (i.e. no other annotation information).
sub length_normalise_matrix
        {
        # Set error messages and accept input parameters
        my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
        my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$matrix_ref)'\n\nwhere".
        "\t\$matrix_ref is a reference to the matrix to be normalised\n\n";
        my @pars = @_ or die $usage;
        foreach my $el (@pars)  {       $el = text::trim($el);  }
        my $matrix_ref = shift @pars or die $usage;      # NB! This will not work if the argument is the number 0 (because it will then be interpreted as false). In that case you need to use 'shift(@pars) or 0', but it may lead to proble$
	my @matrix = @{$matrix_ref};

	# Send the header to an array
	my $headerref = shift(@matrix);
	my @header = @{$headerref};
	my ($length_ind) = grep { $header[$_] eq "length" } 0..$#header; 

	# Loop over @matrix and normalise all values for each gene
	for(my $cc=0; $cc<=$#matrix; $cc++)
		{
		my $length = $matrix[$cc][$length_ind];
		my $norm_factor = 1000/$length;		

		# Loop over sample expression values and normalise them
		for(my $dd=1; $dd<$length_ind; $dd++)
			{
			my $value = $matrix[$cc][$dd];
			$matrix[$cc][$dd] = $value * $norm_factor;
			}
		} 

	# Add the header again
	unshift(@matrix, \@header);
        return(@matrix);
        }

# Strips the Feature_IDs in the first column an expression matrix of the ending feature length information (format "gene_ABC_23_445") and
# converts it into a length column that is added to the matrix
sub length_from_id_to_column
        {
        # Set error messages and accept input parameters
        my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
        my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$infile.csv, \$header_y_n, \$outname)'\n\nwhere".
        "\t\$infile.csv is a gene expression matrix in csv format\n".
	"\t\$header_y_n is an indicator of whether the infile has a header or not ('y' for yes, 'n' for no)\n".
	"\t\$outname is the preferred name of the outfile\n\n";
        my @pars = @_ or die $usage;
        foreach my $el (@pars)  {       $el = text::trim($el);  }
        my $infile = shift @pars or die $usage;      # NB! This will not work if the argument is the number 0 (because it will then be interpreted as false). In that case you need to use 'shift(@pars) or 0', but it may lead to probl$
	my $header_y_n = shift @pars or die $usage;
        my $outname = shift @pars or die $usage;

	# Read the infile
	my @matrix=fileTools::read_table($infile, "csv");	
	my $headerref="";

	if($header_y_n eq "y")	{	$headerref = shift(@matrix);	}
	my $last_ind = scalar(@{$matrix[0]});

	# Loop over @matrix and convert length information to a new column value
	for(my $cc=0; $cc<=$#matrix; $cc++)
		{
		my $id = $matrix[$cc][0];
		my @parts = split("_", $id);
		my $last_part = $#parts;
		my $new_id = join("_", @parts[0..($last_part-2)]);
		my $length = $parts[$last_part] - $parts[($last_part)-1] + 1;

		$matrix[0] = $new_id;
		$matrix[$last_ind] = $length;
		}
	
	if($header_y_n eq "y")
		{
		my @header = @{$headerref};
		push(@header, "length");
		unshift(@matrix, \@header);
		}

	fileTools::write_table(\@matrix, "csv", $outname, "lin");
        }

return(1);

# end functions
