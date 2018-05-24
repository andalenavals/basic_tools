def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='A simple program to find the union of a set of txt files with only one column of ints')
    
    parser.add_argument('--files',  nargs='+', default='', type=str, 
                        help='Files to find the union,  file1 file2 file3...')
    parser.add_argument('--outname', default='', type=str, 
                        help='Name of the output file example zone01.riz')

    args = parser.parse_args()

    return args

def main():
    import numpy as np
 
    args = parse_args()
    setiter= set(np.loadtxt(args.files[0],dtype= 'i4' ))
    for fil in args.files[1:] :
        setaux= set(np.loadtxt(fil,dtype= 'i4' ))
        setiter =  setiter.union(setaux)
    if(args.outname):
        np.savetxt(args.outname, np.array(list(sorted(setiter))), fmt='%i')
    else:
        print(list(sorted(diff)))
         
        
    
if __name__ == "__main__":
    main()

#import itertools as itt
#for files,  i in itt.product(args.files, range(len (args.files)) ):

 
