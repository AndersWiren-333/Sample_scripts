This repository is maintained by Anders Wirén on GitHub
(https://github.com/AndersWiren-333/Sample_scripts)

*Repository structure*
All scripts belonging to this repository are stored in the "scripts" folder, and all modules/
packages used by the scripts are in "modules". The repository has been built so that scripts and
modules can easily find and refer to each other given the folder structure. If you want
to use individual scripts outside the repository folder structure, make sure you manually
set the path to the other scripts and modules (which the scripts may depend on) in the
header of those scripts. The reason for this arrangement (rather than simply putting the scripts
and modules folders in your local path variable) is to avoid accidentally referring to the wrong
modules/scripts when two or more similar repositories (with scripts/modules of equal name) are
present on the same system.

There is also a "maintainance" folder where template scripts for various languages are stored,
along with scripts that aid maintainance and updating of the repository.

The "language_references" folder contains textfiles with syntax and usage examples for common
techniques in the various languages.

*Repository maintenance*
Most Perl scripts have a universal script header which makes it easier to call modules and
other scripts from within them, and which can be updated in all scripts simultaneously by
modifying the header of "template_script.pl" and then running "update_perl_headers.pl".
This script checks for the flag "# perl_script_update" in the header of the relevant scripts.
Similarly, there is a "template_module.pm" and "update_perl_headers.pl" uses this to update
the relevant modules in the same way.

Currently, these update mechanisms are in place for (most) perl scripts, but the aim is to
set them up also for R and SQL scripts.


*Scripts and functions list*
"Scripts_and_functions.xlsx" contains a list of of perl scripts and functions
(subroutines), including descriptions and usage templates that can be copy-pasted
into new scripts as needed. The list is still under construction and incomplete.
When you write a new function, module or script, please update this file!
