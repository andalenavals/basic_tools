## basic_tools
Miscelaneous set of basic tools, useful in data analysis, focused in cosmology. 
All the scripts have the usuals --help and -h flags, if you want more detailed information.
  
  # 1) split.py 
    
  This program divides files in sets of rows, then it generates output files for each set. The maximum difference between the number of rows of all the generated files is 1. It redistribuite the residue over all the sets without changing the order of the original rows.
  
   HOW TO RUN
    
      $ python spliy.py --file=zone.riz --parts=7
  
  
 # 2) fitinfo.py
  
  This program reads fits files. 
  
   HOWTO RUN 
  
   If you want to know the basic informaion inside the file do
    
      $ python fitinfo.py --file=myfile.fits

   Now, if you want to look inside an specific HDU you have to select it with the flag --hdu. Might be you want to know first the header information of the hdu and if it have columns see there names, so the following two additional booleans flags will help.

      $ python fitinfo.py --file=myfile.fits --hdu=1 --header --columns

   Aditional if the hdu is a table you can see all the data of the table,

      $ python fitinfo.py --file=myfile.fits --hdu=1 --alldata 

   or if it have too many columns and you do not want to print all of them, you can select a field
  
      $ python fitinfo.py --file=myfile.fits --hdu=1 --fields 'ra' 'dec' 'mag'
   

# 3) interesct.py

  This program simply, find the intersection between any number of sets. The sets are files with one column of ints. This is useful for instance to compare find the expousres in common of different zones.
  
    $ python intersect --files zone01.riz zone02.riz zone03.riz
  
# 4) maching.py
  
  This program make a match of two files based on their R.A and DEC information. Usually, we want to add one field from one file that is not in the other. The file that have this additional information is call the --fileref and the file of destine just --file. So the outputh after running this script will be, all the information of the --file and the additional --fields included from --fileref. Both files must have fields called 'ra' and 'dec' or any combination of capital or lower letters. 
  
  HOW TO RUN
    
    $ python matching.py --file=psf_cat_232382_1.fits --fileref=COADD_MOF.fits --fields 'CM_T' 'EXT_MOF' --d2d=0.0001 --filename=OUT.fit
  
  After this run the file OUT.fit will have all the information than psf_cat_232382_1.fits and two additional fields 'CM_T' 'EXT_MOF.   
 However, all those rows where the two dimentional separation was higher thatn 1e-4 degree will return None, i.e, they are skipped.
 
  If the --filename is the same than the --file, it will be overwritten.
  
# 5) readcatalog.py
  Each catalogs demands its own version of reading the files, the folling script have 3 examples of how to read catalogs and plot the zone of the sky they cover (Each example are in the function read_alldata, read_alldata2 and read_alldata3).
  
  HOW TO RUN
     
    $ python readcatalog.py --explist=../DESWL/psf/astro/all_zones.riz --fields 'ra' 'dec' --inpath=/home2/dfa/sobreira/alsina/catalogs/y3a1-v23/ --outname=yv23.png

  
  
    
    
