from utils_csn import test_model
import argparse as agp


def main():
    parser = agp.ArgumentParser()
    parser.add_argument('-input','--inputfile',help="specify the inputfile location and name, it has to be in .fasta or .fa format, do not contain letter other than ACGT")
    parser.add_argument('-output','--outputfile',help="the output file location and name")
    args = parser.parse_args()
    
    test_model(args.inputfile,args.outputfile)



if __name__ == '__main__':
    main()
