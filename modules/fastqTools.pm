package fastqTools;

# sub read_fastq($file.fastq)
# sub print_fastq($matrixref)
# sub count_reads_fastq($infile.fastq)
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

# Reads a given fastq file into a matrix (one row = one sequence, four columns per row)
sub read_fastq
        {
	my $subname = "read_fastq";
        my $usage="\nSyntax error for sub ${subname}. Correct usage: '\${packname}::\${subname}(\$infile.fastq)'\n\nwhere".
        "\t\$infile.fastq is the input file, in fastq format\n\n";
        my @pars = @_ or die $usage;
        foreach my $el (@pars)  {       $el = text::trim($el);  }
        my $infile = shift @pars;
	
	# Open file
	open(my $in, "<", $infile) or die "sub $subname couldn't open infile $infile\n";

        #if(!open($in, "<", $infile))   {       die "sub $subname couldn't open infile $infile\n";    }
        my @matrix=();
        my @temp=();
        my $counter=0;
        while(my $line = <$in>)
                {
                $counter++;
                $line = text::trim($line);
                push(@temp, $line);
                if($counter==4)                         # If this is the last line of a sequence record...
                        {
                        push(@matrix, [@temp]);         # Add the four lines of the record as an element (a line of) to @one matrix
                        @temp=();                                       # Free the @temp array
                        $counter=0;                             # Reset the counter
                        }
                }
	close($in);
        return(@matrix);
        }

# Prints a given fastq matrix (one row = one sequence, four columns per row) to a file with the specified name
sub print_fastq
	{
	# Set error messages and accept input parameters
        my $subname = "print_fastq";
        my $usage="\nSyntax error for sub ${subname}. Correct usage: '\${packname}::\${subname}(\$matrixref, \$outname)'\n\nwhere".
        "\t\$matrixref is a reference to a matrix holding fastq reads (one row = one read, each row has four columns)\n".
	"\t\$outname is the deired name of the fastq outfile\n\n";
        my @pars = @_ or die $usage;
        foreach my $el (@pars)  {       $el = text::trim($el);  }
	my $matrix_ref = shift(@pars) or die $usage;
        my $outname = shift(@pars) or die $usage;

	# Start doing things
        my @matrix = @{$matrix_ref};
	open(my $out, ">>", $outname) or die "sub $subname couldn't create outfile $outname\n";

        for(my $c=0; $c<=$#matrix; $c++)
                {
                my @arr=@{$matrix[$c]};
                foreach my $el (@arr)
			{
			$el = text::trim($el);
			print($out "$el\n");
			}
                }
        close($out);
	}

# Counts the number of sequence reads in a fastq file and returns the number.
sub count_reads_fastq
	{
	# Set error messages and accept input parameters
        my ($calling_script, $calling_line, $subname) = (caller(0))[1,2,3];
        my $usage="\nUsage error for subroutine '${subname}' (called by script '${calling_script}', line ${calling_line}). Correct usage: '${subname}(\$infile.fastq)'\n\nwhere".
        "\t\$infile.fastq is the file to have its reads counted\n\n";
        my @pars = @_ or die $usage;
        foreach my $el (@pars)  {       $el = text::trim($el);  }
        my $infile = shift @pars or die $usage;
	open(my $in, "<", $infile) or die "Sub $subname could't open infile $infile\n";

	my $lines=0;
	while(my $line = <$in>)
		{
		$lines++;
		}
	close($in);
	my $seqs = $lines/4;
	if(int($seqs) != $seqs)
		{
		die "Subroutine $subname finds that the number of lines in fastq file $infile is not evenly divisible with 4. This indicates there is something wrong with the formatting of the file\n";
		}
        return($seqs);
	}

return(1);

# end functions
