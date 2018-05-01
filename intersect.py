def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='A simple program to find the intesection of a set of txt files with only one column of ints')
    
    parser.add_argument('--files',  nargs='+', default='', type=str, 
                        help='Files to find the intersection')
    parser.add_argument('--outname', default='', type=str, 
                        help='Name of the output image example zone01.png')

    args = parser.parse_args()

    return args

def main():
    import numpy as np
 
    args = parse_args()
    inter= set(np.loadtxt(args.files[0],dtype= 'i4' ))
    for files in args.files[1:]:
        print(files)
        setaux= set(np.loadtxt(files,dtype= 'i4' ))
        inter =  inter.intersection(setaux)
    print(list(sorted(inter)) )
    if(args.outname):
        np.savetxt(args.outname, np.array(list(sorted(inter))), fmt='%i')
         
        
    
if __name__ == "__main__":
    main()

#import itertools as itt
#for files,  i in itt.product(args.files, range(len (args.files)) ):

 
