def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='A simple program to find the intesection of a set of txt files with only one column of ints')
    
    parser.add_argument('--files',  nargs='+', default='', type=str, 
                        help='Files to find the intersection')

    args = parser.parse_args()

    return args

def main():
    import numpy as np
 
    args = parse_args()
    inter= set(np.loadtxt(args.files[0],dtype= 'i4' ))
    for files in args.files:
        setaux= set(np.loadtxt(files,dtype= 'i4' )).intersection()
        inter =  inter.intersection(setaux)
    print(inter)
         
        
    
if __name__ == "__main__":
    main()

#import itertools as itt
#for files,  i in itt.product(args.files, range(len (args.files)) ):

 
