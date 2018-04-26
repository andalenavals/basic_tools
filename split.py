def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Split files program')
    
    parser.add_argument('--file', default='',
                        help='File to split')
    parser.add_argument('--parts', default='',
                        help='Number of parts')
    args = parser.parse_args()

    return args

def testinputfiles(args):
    if args.file != '':
        print ('\n File ',args.file, 'sucessfully read. OK! \n' )
    else:
        print ('No input file1')


def main():
    import fitsio
    from astropy.io import fits
    import numpy as np

    args = parse_args()
    testinputfiles(args)
    
    dt =  np.dtype([('exposure','>i4' )])
    data = np.loadtxt(args.file,  dtype=dt)
    parts =  int(args.parts) 
    newlen = len(data) // parts
    res = len(data) % parts

    if(res!=0):
        lowlimit = 0
        uplimit =  (newlen + 1)
        for k in range (res):
            print(k ,  lowlimit ,  uplimit)
            name =  args.file[:- 4] + '_' + str(k + 1) + args.file[- 4: ]
            np.savetxt(name, data[lowlimit : uplimit], fmt='%i')
            lowlimit = uplimit
            uplimit =  lowlimit + (newlen + 1)
        lowlimit2 = lowlimit
        uplimit2 = lowlimit2 + newlen 
        for k in range (res, parts):
            print(k,  lowlimit2,  uplimit2)
            name =  args.file[:- 4] + '_' + str(k + 1) + args.file[- 4: ]
            np.savetxt(name, data[lowlimit2: uplimit2], fmt='%i')
            lowlimit2 = uplimit2
            uplimit2 =  lowlimit2 + newlen
    else:
        for k in range (parts):
            name =  args.file[:- 4] + '_' + str(k + 1) + args.file[- 4: ]
            np.savetxt(name, data[k*(newlen) : (k + 1) * newlen ], fmt='%i')
    

    

if __name__ == "__main__":
    main()
    




    #data = fitsio.read(args.file)
    #print(data)
    #wsc_range(data)
    
    
 
