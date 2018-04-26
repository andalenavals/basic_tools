# basic_tools
Miscelaneous set of basic toolS, useful in data analysis, focused in cosmology.
All the scripts have the usuals --help and -h flags, if you want more detailed information.

1) split.py
This program divides files in sets of rows, then it generates output files for each set. The maximum difference between the number
of rows of all the generated files is 1. It redistribuite the residue over all the sets without changing the order of the original
rows.
  HOW TO RUN
  $ python spliy.py --file=zone.riz --parts=7
  
