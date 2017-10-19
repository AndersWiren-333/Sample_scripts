This repository is maintained by Anders Wirén on GitHub
(https://github.com/AndersWiren-333/Sample_scripts)

*Repository structure*
All scripts belonging to this repository are stored in the "scripts" folder, and all modules
used by the scripts are in "modules". The repository has been built so that scripts and
modules can easily find and refer to each other given this folder structure. If you want
to use individual scripts outside the repository folder structure, make sure you manually
set the path to the other scripts and modules (which the scripts may depend on) in the
header of those scripts.

*Repository maintenance*
Scripts with a version number have the universal script header (for the respective language,
e.g. Perl or R) which makes it easier to call modules and other scripts from within them,
and which can be updated in all scripts simultaneously by modifying the header of
"template_script.pl" and then running "update_universal_script_header.pl". This script uses
"scripts/script_update_list.txt" to determine which scripts should be updated and which
should not (to exclude a script from updating, either leave it out of the list or put a "#"
before its name in the list). Similarly, there is a "template_module.pm",
"update_universal_module_header.pl" and "module_update_list.txt" for updating module
headers.

Currently, these update mechanisms are in place for (most) perl scripts, but the aim is to
set them up also for R scripts and other languages that may be used in the future.


*Modules and functions list*
"Modules_and_functions.xlsx" contains a list of of perl modules and their functions
(subroutines), including descriptions and usage templates that can be copy-pasted
into new scripts as needed. A similar list for scripts is under construction.
When you write a new function, module or script, please update this file!
