package fastaTools;

# sub fasta_to_matrix($fasta_file)
# sub fasta_to_matrix2($fasta_file)
# sub divide_fasta($fasta_file, $num_seqs)
# sub seqs_from_reference(\@feat_ids, $ref_sequence)
# sub count_nr_fasta_reads($fasta_file)
# sub count_reads_fasta($fasta_file)
# sub count_reads_mnr_fasta($fasta_file)
# sub deredundise_nr_fasta($file, $outname)
# sub get_seqs_from_fasta($ID_list.txt, $fasta_file.fa, $outname)
# sub get_seq_lengths_from_fasta($infile.fa);
# sub remove_alignment_target_from_fasta($infile.fa, $target.fa, $mismatch, $outname)
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

# Reads a fasta file into a matrix
sub fasta_to_matrix
	{
	# Set error messages and accept input parameters
	my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
	my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$fasta_file)'\n\nwhere".
	"\t\$fasta_file is the file to read into a matrix\n\n";
	my @pars = @_ or die $usage;
	foreach my $el (@pars)  {       $el = text::trim($el);  }
	my $infile = shift @pars or die $usage;

	my $fa = Bio::SeqIO->new('-format' => 'fasta', '-file' =>  "$infile");
	my @matrix=();

	while(my $seq = $fa->next_seq())
		{
        my $sequence = $seq->seq();
		$sequence = uc(text::trim($sequence));
        my $id = $seq->display_id();
		$id =~ s/>//;
		$id = text::trim($id);
		push(@matrix, [$id, $sequence]);
		}

	return(@matrix);
	}

# Reads a fasta file into a matrix
sub fasta_to_matrix2
	{
        # Set error messages and accept input parameters
        my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
        my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$fasta_file)'\n\nwhere".
        "\t\$fasta_file is the file to read into a matrix\n\n";
        my @pars = @_ or die $usage;
        foreach my $el (@pars)  {       $el = text::trim($el);  }
        my $infile = shift @pars or die $usage;

        if(!open(IN, $infile))  {       die "Subroutine $subname (called by script '${calling_script}', line ${calling_line}) couldn't open infile $infile\n"; }
        my @arr=();

        # First add each line as an element to an array
        while(<IN>)
			{
			my $line = text::trim($_);
			push(@arr, $line);
			}
        close(IN);

        # If the first line isn't an id line, exit with an error message
        if($arr[0] !~ m/^>/)    {       die "The fasta file $infile used in sub fasta_to_matrix seems to have the wrong format (it doesn't start with a proper id line). Correct the file and try again.\n";    }
        if($arr[$#arr] =~ m/^>/)        {       die "The fasta file $infile used in sub fasta_to_matrix seems to have the wrong format (it ends with an id line). Correct the file and try again.\n";   }

        # Then add all lines in the correct places in a matrix (i.e. id lines in first column and sequences in the second)
        my @matrix=();
        my $seqno=0;
        my $tempstring="";

        for(my $cc=0; $cc<=$#arr; $cc++)
			{
			my $line = $arr[$cc];

			# If this is an id line
			if($line =~ /^>/)
				{
				# If it is the first id line, just add it to the matrix
				$line =~ s/>//;         # Remove the angle bracket from the ID line
				if($cc==0)      {       $matrix[$seqno][0]=$line;       }

				# If it is any other id line, finish adding the preceding sequence to the matrix, then add the new id line
				else
					{
					$matrix[$seqno][1] = $tempstring;
					$tempstring="";
					$seqno++;
					$matrix[$seqno][0]=$line;
					}
				}

			# Else, if it is a sequence line
			else
				{
				# If it is the last line in the file
				if($cc==$#arr)
					{
					$tempstring = "$tempstring"."$line";
					$matrix[$seqno][1]=$tempstring;
					$tempstring="";
					}

				# If it is any other sequence line
				else    {       $tempstring = "$tempstring"."$line";    }
				}
			} # Loop over elements in array ends here

        return(@matrix);
	}


# DIVIDES A FASTA FILE INTO SMALLER PIECES. The pieces are defined by a given number of sequences per piece
sub divide_fasta
	{
	# Declare variables and filehandles
	my $usage="\nSyntax error for sub divide_fasta. Correct usage: 'fastaTools::divide_fasta_file(\$fasta_file, \$sequences_per_piece)'\n";
	my $infile = $_[0] or die $usage;
	my $num_seqs= $_[1] or die $usage;
	if($num_seqs%2 != 0)	{	die "sub divide_fasta can't handle an uneven number of sequences! Try again.\n";	}
	if(!open(COUNT, $infile))	{	die "sub divide_fasta couldn't open $infile for counting lines\n";	}

	# First count the number of lines in the file
	my $num_lines=0;
	my $first_line="";
	my $last_line="";
	my $seqs_in_file = 0;

	while(<COUNT>)
		{
		$num_lines++;
		my $line = text::trim($_);
		if($line =~ m/^>/)	{	$seqs_in_file ++;	}
		if($num_lines==1)	{	$first_line = $line;	}
		$last_line = $line;
		}
	close(COUNT);

	# Check that the first line is a proper id line and the last line isn't
	if($first_line !~ m/^>/)	{	die "The fasta file $infile used in sub divide_fasta seems to have the wrong format (it doesn't start with a proper id line). Correct the file and try again.\n";	}
	if($last_line =~ m/^>/)	{	die "The fasta file $infile used in sub divide_fasta seems to have the wrong format (it ends with something that looks like an id line). Correct the file and try again.\n";	}

	# Read infile and subdivide it as you go along
	if(!open(IN, $infile))       {       die "sub divide_fasta couldn't open $infile for subdividing\n";      }
	my $seqno=0;
	my $lineno=0;
	my $part=1;
	my $limit=$num_seqs;
	my $tempstring = "";
	#if(!open(OUT, ">>${infile}_part_${part}.fa"))	{	die "Couldn't create outfile\n";	}
	if(!open(OUT, ">>${infile}_${part}.fa"))  {       die "Couldn't create outfile\n";        }

	while(<IN>)
		{
		$lineno++;
		my $line=text::trim($_);		

		# If id line
		if($line =~ m/^>/)
			{
			# If first id line in file/part, just print the line to the outfile
			if($seqno==0)
				{
				$seqno++;
				print(OUT "$line\n");
				}

			# If any other id line
			elsif($seqno<$limit)
				{
				# Finish printing out the previous sequence line
				print(OUT "$tempstring\n");
				$tempstring="";
				$seqno++;		# Indicate that we start a new sequence
				print(OUT "$line\n");	# Print line to outfile
				}

			elsif($seqno==$limit)
				{
				print(OUT "$tempstring\n");
                                $tempstring="";
				$seqno=1;
				close(OUT);
				$part++;
				if(!open(OUT, ">>${infile}_${part}.fa"))    {       die "Couldn't create outfile\n";        }
				print(OUT "$line\n");
				}
			}

		# If sequence line
		else
			{
			# If this is not the last line of the infile
			if($lineno<$num_lines)	{	$tempstring = "$tempstring"."$line";	}

			# If it is the last line of the infile
			elsif($lineno==$num_lines)
				{
				$tempstring = "$tempstring"."$line";
				print(OUT "$tempstring\n");
				close(OUT);
				$tempstring="";
				}
			} # End of else statement
		} # End of second while loop
	close(IN);
	return($part);
	} # End of sub


# Gets the reference sequence for given feature IDs from a given reference genome/transcriptome and prints them to a fasta file
sub seqs_from_reference
	{
	my $usage="\nSyntax error for sub seqs_from_reference. Correct usage: 'fastaTools::seqs_from_reference(\@feat_ids, \$ref_sequence.fa, \$outfile_name)'\n\nwhere".
	"\t\@feat_ids is an array with feature IDs in the format NC_1234.5_678_1011\n".
	"\t\$ref_sequence.fa is a fasta file containing the refernce genome/transcriptome/contigs\n".
	"\t\$outfile_name is what you want the outfile to be called. As yet, this argument has to be given even if you want the results as a matrix\n\n";
	my @pars=@_;
	my $idsref = text::trim($pars[0]);
	my $ref_genome = text::trim($pars[1]);
	my $outname = text::trim($pars[2]);
	my @ids = @{$idsref};
	if(!open(OUT, ">>$outname"))	{	die "sub seqs_from_reference couldn't create outfile gene_reference_seqs.fa\n";	}	

	# Read reference genome/transcriptome/sequence into matrix;
	my @reference = fasta_to_matrix($ref_genome);
	
	my @matrix=();
	# Loop over feature ids
	for(my $cc=0; $cc<=$#ids; $cc++)
		{
		my $id = $ids[$cc];
		my @id_parts = split("_", $id);	
		my $seqid = "$id_parts[0]"."_"."$id_parts[1]";
		my $start = $id_parts[2];
		my $stop = $id_parts[3];
		my $length = ($stop-$start)+1;
		
		# Loop over fasta reference sequences
		REFS: for(my $dd=0; $dd<=$#reference; $dd++)
			{
			my @refid_parts = split(" ", $reference[$dd][0]);
			my $refid = text::trim($refid_parts[0]);
			$refid =~ s/>//;	# Remove the angle bracket from the fasta ID line , if it is there
			if($refid eq $seqid)
				{
				my $refseq = $reference[$dd][1];
				my $new_seq = substr($refseq, $start, $length);
				print(OUT ">$id\n$new_seq\n");
				push(@matrix, [$id, $new_seq]);
				last REFS;
				}
			}
		}
	close(OUT);
	return(@matrix);
	}

# Counts the number of total and unique (= non-redundant) reads in a fasta file in non-redundant format
sub count_nr_fasta_reads
	{
	# Set error messages and accept input parameters
        my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
        my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$infile.fa)'\n\nwhere".
        "\t\$infile.fa is the file (in non-redundant fasta format) which should have its reads counted\n\n";
	my @pars = @_ or die $usage;
        foreach my $el (@pars)  {       $el = text::trim($el);  }
        my $infile = shift @pars or die $usage;

	open(my $in, "<", $infile) or die "Subroutine '$subname' (called by script '${calling_script}', line ${calling_line}) couldn't open infile $infile\n";

	my $total=0;
	my $unique=0;

	# Loop over reads in file
	while(my $line = <$in>)
        	{
        	my $line = text::trim($line);
        	if($line =~ m/>/)
                	{
                	$unique++;
                	my @arr = split("-", $line);
			foreach my $el (@arr)	{	$el = text::trim($el);	}
                	my $abund = $arr[1];
                	$total = $total+$abund;
                	}
        	}
	close($in);
	my @stats = ($total, $unique);
	return(@stats);
	}

# Counts the number of sequence reads in (any) fasta file and returns that number
sub count_reads_fasta
	{
	# Set error messages and accept input parameters
        my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
        my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$fasta_file.fa)'\n\nwhere".
        "\t\$fasta_file.fa is the file to have its reads counted\n\n";
        my @pars = @_ or die $usage;
        foreach my $el (@pars)  {       $el = text::trim($el);  }
        my $infile = shift @pars or die $usage;

	open(my $in, "<", $infile) or die "Subroutine $subname (called by script '${calling_script}', line ${calling_line}) couldn't open infile $infile\n";
	my $id_lines=0;

	while(my $line = <$in>)
		{
		$line = text::trim($line);
		if($line =~ m/^>/)	{	$id_lines++;	}
		}
	close($in);
        return($id_lines);
	}

# Counts the number of total and redundant sequence reads in a fasta file in multi_redudnant format, i.e.
# >sequence_abundunce1_abundance2_abundance3...
# sequence (where each abundance in the id line is the abundance of that read in a specific sample/library)
sub count_reads_mnr_fasta
        {
        # Set error messages and accept input parameters
        my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
        my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$infile.fa)'\n\nwhere".
        "\t\$infile.fa is the file whose reads should be counted\n\n";
        my @pars = @_ or die $usage;
        foreach my $el (@pars)  {       $el = text::trim($el);  }
        my $infile = shift @pars or die $usage;

	# Read file into matrix (two columns: id and sequence)
	my @matrix=fastaTools::fasta_to_matrix($infile);
	my $total=0;
	my $unique=0;

	for(my $cc=0; $cc<=$#matrix; $cc++)
		{
		$unique++;
		my @id_parts = split("_", $matrix[$cc][0]);	
		shift(@id_parts);
		for(my $dd=0; $dd<=$#id_parts; $dd++)
			{
			#$id_parts[$dd] = text::trim($id_parts[$dd]);
			$total = $total+$id_parts[$dd];
			}
		}
	my @stats = ($total, $unique);
        return(@stats);
        }

# Deredundises a fasta file in non-redundant format (i.e. '>sequnce-abudance_linebreak_sequence')
sub deredundise_nr_fasta
	{
        # Set error messages and accept input parameters
        my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
        my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$infile.fa, \$outname)'\n\nwhere".
        "\t\$infile.fa is the file to be deredundised\n".
	"\t\$outname is the desired name of the resulting file\n\n";

        my @pars = @_ or die $usage;
        foreach my $el (@pars)  {       $el = text::trim($el);  }
        my $infile = shift @pars or die $usage;      # NB! This will not work if the argument is the number 0 (because it will then be interpreted as false). In that case you need to use 'shift(@pars) or 0', but it may lead to proble$
	my $outname = shift @pars or die $usage;

	# Read the file into a matrix
	my @matrix=fastaTools::fasta_to_matrix($infile);

	# Sort on matrix on sequence
	@matrix = sort { $a->[1] cmp $b->[1] } @matrix;

	# Loop over matrix and de-redundise
	my @new_matrix=();
	my $compseq=$matrix[0][1];
	my ($compsekv, $comp_abund) = split("-", $matrix[0][0]);

	for(my $cc=1; $cc<=$#matrix; $cc++)
		{
		my $seq = $matrix[$cc][1];
		my ($sekv, $abund) = split("-", $matrix[$cc][0]);
	
		# If we are still on the same sequence
		if($seq eq $compseq)
			{
			$comp_abund += $abund;

			# If this is the last line of the matrix
			if($cc==$#matrix)
				{
				my $new_id = ">".$compseq."-".$comp_abund;
				push(@new_matrix, [($new_id, $compseq)]);
				}
			}

		# If we are on a new sequence
		else
			{
			my $new_id = ">".$compseq."-".$comp_abund;
			push(@new_matrix, [($new_id, $compseq)]);

			$compseq = $seq;
			$comp_abund = $abund;

			# If this is the last line of the matrix
			if($cc==$#matrix)
				{
				$new_id = ">".$compseq."-".$comp_abund;
                        	push(@new_matrix, [($new_id, $compseq)]);
				}
			}
		} # Loop over original matrix ends

	# Print out the new matrix
	fileTools::write_table_as_list(\@new_matrix, $outname);
	}


# Given a list of sequence IDs (as a textfile with one ID on each line), retrieves the sequences of those IDs from a fasta file
# and prints the selected IDs and sequences to a new fasta file of the specified name.
sub get_seqs_from_fasta
        {
        # Set error messages and accept input parameters
        my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
        my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$ID_list.txt, \$fasta_file.fa, \$outname)'\n\nwhere".
	"\t\$ID_list.txt is a list of sequence IDs\n".
        "\t\$fasta_file.fa is the fasta files to retrieve sequence from\n".
        "\t\$outname is the desired name of the resulting fasta file\n\n";

        my @pars = @_ or die $usage;
        foreach my $el (@pars)  {       $el = text::trim($el);  }
        my $id_list = shift @pars or die $usage;      # NB! This will not work if the argument is the number 0 (because it will then be interpreted as false). In that case you need to use 'shift(@pars) or 0', but it may lead to proble$
	my $fasta_file = shift @pars or die $usage;
        my $outname = shift @pars or die $usage;

	# Read the ID list into an array
	my @ids=misc::read_list($id_list);
	@ids = sort { $a cmp $b } @ids;

        # Read the fasta file into a matrix
        my @matrix=fastaTools::fasta_to_matrix($fasta_file);

        # Sort on matrix on ID
        @matrix = sort { $a->[0] cmp $b->[0] } @matrix;

	# Pick out the sequences of the listed IDs from the fasta file and write the selected IDs and sequences
	# to a new fasta file
	my @sel_matrix=();

	# Loop over selected IDs
	IDS: for(my $cc=0; $cc<=$#ids; $cc++)
		{
		# Loop over sequences in @matrix
		for(my $dd=0; $dd<=$#matrix; $dd++)
			{
			if($ids[$cc] eq $matrix[$dd][0])
				{
				my $new_id = ">"."$ids[$cc]";
				push(@sel_matrix, [($new_id, $matrix[$dd][1])]);
				next IDS;
				}
			}
		} # Loop over @ids ends

	# Print the new @sel_matrix to a fasta outfile
	fileTools::write_table_as_list(\@sel_matrix, $outname);
        }

# Takes a fasta file as input, computes the length of each sequence in it, and returns a matrix where each row has two columns: a sequence identifier
# and the corresponding length.
sub get_seq_lengths_from_fasta
        {
        # Set error messages and accept input parameters
        my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
        my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$infile)'\n\nwhere".
        "\t\$infile.fa is an input file in fasta format\n\n";
        my @pars = @_ or die $usage;
        foreach my $el (@pars)  {       $el = text::trim($el);  }
        my $infile = shift @pars or die $usage;      # NB! This will not work if the argument is the number 0 (because it will then be interpreted as false). In that case you need to use 'shift(@pars) or 0', but it may lead to proble$

	my @outmatrix=();
	push(@outmatrix, [("Feature_ID", "length")]);

	# Read the fasta file
	my @fasta_matrix=fastaTools::fasta_to_matrix($infile);
	
	# Loop over @fasta_matrix and compute lengths
	for(my $cc=0; $cc<=$#fasta_matrix; $cc++)
		{
		my $id = $fasta_matrix[$cc][0];
		my $len = length($fasta_matrix[$cc][1]);
		push(@outmatrix, [($id, $len)]);
		}

        return(@outmatrix);
        }

# Aligns the sequence reads in a fasta file (the infile) to those in a different fasta file (the target) using Patman (Prüfer et al 2008), with a given number of mismatches allowed between
# infile sequence and target sequence. Those reads from the infile that align to something in the target are removed from the infile. The original infile is preserved, and a new infile with
# the preferred name is produced. The infile can be of either redundant ("r", i.e. '>id (linebreak) sequence'), non-redundant ("nr", i.e. '>sequence-abundance (linebreak) sequence')
# or multi-non-redundant ("mnr", '>sequence_abundance1_abundance2_abundance3... (linebreak) sequence') format. Also specify the maximum the number of sequences that the infile can contain before
# it gets split into smaller pieces before alignment, to save memory, and whether the infile, the outfile or both should be compressed after use.
sub remove_alignment_target_from_fasta
        {
        # Set error messages and accept input parameters
        my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
        my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$infile.fa, \$target.fa, \$mismatch, \$outname,".
	" \$seqs_per_infile, \$infile_format, \$zip_in_out_all)'\n\nwhere".
        "\t\$infile.fa is the fasta file that something should be removed from\n".
        "\t\$target.fa is the fasta file to align the infile to\n".
	"\t\$mismatch is the number of mismatches to allow in the alignment between infile and target\n".
        "\t\$outname is the preferred name of the output file\n".
	"\t\$seqs_per_infile is the number of sequences that the infile can contain before it gets split into smaller pieces before alignment (the pieces will be put together again in the end)\n".
	"\t\$infile format can be either 'r', 'nr' or 'mnr' (redundant, non-redundant, multi-non-redundant), where\n".
	"\t\t'r' is 'redundant' i.e. '>id (linebreak) sequence'\n".
	"\t\t'nr' is 'non-redundant' i.e. '>sequence-abundance (linebreak) sequence'\n".
	"\t\t'mnr' is 'multi-non-redundant' i.e. '>sequence_abundance1_abundance2_abundance3... (linebreak) sequence'\n".
	"\t\$zip_in_out_all is an indicator of whether either input or output files or both should be compressed after the script has finished (options are 'in', 'out' or 'all'. Default is no compression)\n";

        my @pars = @_ or die $usage;
        foreach my $el (@pars)  {       $el = text::trim($el);  }
        my $infile = shift @pars or die $usage;      # NB! This will not work if the argument is the number 0 (because it will then be interpreted as false). In that case you need to use 'shift(@pars) or 0', but it may lead to proble$
	my $target = shift @pars or die $usage;
	my $mm = shift @pars or die $usage;
	my $outname = shift @pars or die $usage;
	my $seqs_per_infile = shift @pars or die $usage;
	my $infile_format = shift @pars or die $usage;

	# If $zip is not defined, set it to "none"
	my $zip="none";
	if($pars[0])	{	$zip = $pars[0];	}

	# Uncompress infile and set file names
	my ($file_gz, $file)=misc::check_compressed($infile);
	my ($dir, $file2, $basename, $suffix)=text::parse_path($file);
	my $patfile = $basename . ".pat";
	
	# Align
	my ($countT, $countU, $complex, $cop_per_read, $countM)=misc::patman_align($infile, $patfile, $target, $seqs_per_infile, $mm, $infile_format);

	# Read the patman outfile
	my @pat_matrix=fileTools::read_table($patfile, "tsv");

	# Pick out column 2
	my @match_reads=matrixTools::get_matrix_columns_as_rows(\@pat_matrix, 1);

	# De-redundize @match_reads
	my @nr_match_reads=misc::unique_list(\@match_reads, "alph");

	# Read the fasta infile
	my @fasta_matrix=fastaTools::fasta_to_matrix($file);

	open(my $out, ">", $outname) or die "Subroutine $subname (called by script '${calling_script}', line ${calling_line}) couldn't create outfile $outname\n";

	# Loop over @fasta_matrix
	for(my $cc=0; $cc<=$#fasta_matrix; $cc++)
		{
		my $keep="y";

		# If there are matching reads, loop over @nr_match_reads
		if(defined($nr_match_reads[0]))
			{
			for(my $dd=0; $dd<=$#nr_match_reads; $dd++)
				{
				if($fasta_matrix[$cc][0] eq $nr_match_reads[$dd])	{	$keep = "n";	}
				}
			unless($keep eq "n")	{	print($out ">$fasta_matrix[$cc][0]\n$fasta_matrix[$cc][1]\n");	}
			}
		else	{	print($out ">$fasta_matrix[$cc][0]\n$fasta_matrix[$cc][1]\n");	}
		}

	unlink($patfile);

	if($zip eq "in")	{	gzip "$infile" => "${infile}.gz"; unlink($infile);	}
	elsif($zip eq "out")	{	gzip "$infile" => "${infile}.gz"; unlink($infile);	}
	elsif($zip eq "all")
		{
		gzip "$infile" => "${infile}.gz"; unlink($infile);
		gzip "$infile" => "${infile}.gz"; unlink($infile);
		}

	close($out);
	#return($something);
	}

# Description of what the function does and how to use it
sub fasta_from_patman
	{
	# Set error messages
	my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
	my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$patfile, \$outname)'\n\nwhere".
	"\t\$patfile is the name of the patman file to be converted to fasta\n".
	"\t\$outname is the preferred name of the fasta outfile\n\n";
	
	# Accept input parameters
	my @pars = @_ or die $usage;
	foreach my $el (@pars)  {       $el = text::trim($el);  }
	my $infile = shift @pars or die $usage;	
	my $outname = shift @pars or die $usage;

	# Uncompress infile if necessary
	my ($patfile_gz, $patfile)=misc::check_compressed($infile);
	
	# Read infile into matrix (patman outfiles are in tsv format)
	my @matrix=fileTools::read_table($patfile, "tsv");

	# Get sequences
	my @seqs=matrixTools::get_matrix_columns_as_rows(\@matrix, 2);
	
	# Remove any duplicates from @seqs
	my @useqs=misc::unique_list(\@seqs, "alph");

	# Loop over sequences and print them to an outfile
	open(my $out, ">", $outname) or die "Subroutine $subname (called by script '${calling_script}', line ${calling_line}) couldn't create outfile $outname\n";
	
	for(my $cc=0; $cc<=$#useqs; $cc++)
		{
		my $id = $useqs[$cc];
		my ($seq, $abund) = split("-", $id);
		print $out "$id\n$seq\n";
		}
		
	close($out);
	
	my $num=scalar(@useqs);
	print "\n$num\n";
	}	
	
return(1);

# end functions
