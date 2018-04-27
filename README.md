# basic_tools
Miscelaneous set of basic tools, useful in data analysis, focused in cosmology.
All the scripts have the usuals --help and -h flags, if you want more detailed information.

1) split.py
This program divides files in sets of rows, then it generates output files for each set. The maximum difference between the number
of rows of all the generated files is 1. It redistribuite the residue over all the sets without changing the order of the original
rows.
  HOW TO RUN
    
    $ python spliy.py --file=zone.riz --parts=7
  
  
2) fitinfo.py
This program reads fits files. 
  HOWTO RUN 
  If you want to know the basic informaion inside the file do
    
    $ python fitinfo.py --file=myfile.fits

  Now, if you want to look inside an specific HDU you have to select it with the flag --hdu. Might be you want to know first the 
  header information of the hdu and if it have columns see there names, so the following two additional booleans flags will help.

    $ python fitinfo.py --file=myfile.fits --hdu=1 --header --columns

  Aditional if the hdu is a table you can see all the data of the table,

    $ python fitinfo.py --file=myfile.fits --hdu=1 --alldata 

  or if it have too many columns and you do not want to print all of them, you can select a field
  
    $ python fitinfo.py --file=myfile.fits --hdu=1 --fields 'ra' 'dec' 'mag'
   
   

  
