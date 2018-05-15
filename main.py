# Libraries
import argparse

# Modules
import average

# Defined constants
ESBL_AB_RESISTENCE_LIST = ["cefotaxim", "ceftazidim", "ceftriaxon"]
CULTURE_SIZE_CUTOFF = 5
AB_CULTURE_COUNT_CUTOFF = 20


def main(args):
    if args.testfile:
        filename = "data/sample.csv"
    elif args.file:
        filename = args.file
    else:
        print("ERROR: Invalid or no filename")
        return
    if args.average:
        average.run(filename, args.culture_cutoff, args.ab_count_cutoff, ESBL_AB_RESISTENCE_LIST)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-file', type=str, help='File to use for analysis')
    parser.add_argument('-testfile', action='store_true', help='Use data/sample.csv')
    parser.add_argument('-average', action='store_true', help='Print results from baseline analysis method')
    parser.add_argument('-culture_cutoff', type=int, help='Used by average: Determines the minimum amount of cultures per patient', default=5)
    parser.add_argument('-ab_count_cutoff', type=int, help='Used by average: Determines the minimum count of AB culture occurences', default=20)
    main(parser.parse_args())