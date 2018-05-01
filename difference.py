def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='A simple program to find the difference of a set of txt files with only one column of ints')
    
    parser.add_argument('--files',  nargs='+', default='', type=str, 
                        help='Files to find the intersection, remember that substraction is not conmutative, so the order of the input matters file1-file2-file3...')
    parser.add_argument('--outname', default='', type=str, 
                        help='Name of the output image example zone01.png')

    args = parser.parse_args()

    return args

def main():
    import numpy as np
 
    args = parse_args()
    diff= set(np.loadtxt(args.files[0],dtype= 'i4' ))
    for fil in args.files[1:] :
        setaux= set(np.loadtxt(fil,dtype= 'i4' ))
        diff =  diff.difference(setaux)
    print(list(sorted(diff)))
    if(args.outname):
        np.savetxt(args.outname, np.array(list(sorted(diff))), fmt='%i')
         
        
    
if __name__ == "__main__":
    main()

#import itertools as itt
#for files,  i in itt.product(args.files, range(len (args.files)) ):

 
