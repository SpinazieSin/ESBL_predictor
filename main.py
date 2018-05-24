# Libraries
import argparse

# Modules
import average
import SVM
import decision_tree
import perceptron

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
    testmode = args.testmode
    if args.average:
        method = average.Average(filename, args.culture_cutoff, args.ab_count_cutoff, ESBL_AB_RESISTENCE_LIST)
    elif args.svm:
        SVM.run(filename, args.culture_cutoff, args.ab_count_cutoff, ESBL_AB_RESISTENCE_LIST)
    elif args.tree:
        decision_tree.run(filename, args.culture_cutoff, args.ab_count_cutoff, ESBL_AB_RESISTENCE_LIST, testmode=testmode)
    elif args.perceptron:
        method = perceptron.Perceptron(filename, args.culture_cutoff, args.ab_count_cutoff, ESBL_AB_RESISTENCE_LIST, testmode=testmode)
    method.run()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-file', type=str, help='File to use for analysis')
    parser.add_argument('-testfile', action='store_true', help='Use data/sample.csv')
    parser.add_argument('-average', action='store_true', help='Train and print results from baseline analysis method')
    parser.add_argument('-svm', action='store_true', help='Train and print results from SVM method')
    parser.add_argument('-tree', action='store_true', help='Train and print results from decision tree method')
    parser.add_argument('-perceptron', action='store_true', help='Train and print results from multilayer perceptron neural network')
    parser.add_argument('-culture_cutoff', type=int, help='Used by average: Determines the minimum amount of cultures per patient', default=5)
    parser.add_argument('-ab_count_cutoff', type=int, help='Used by average: Determines the minimum count of AB culture occurences', default=20)
    parser.add_argument('-testmode', type=str, help='Either test_patient or cross_validation', default='cross_validation')
    main(parser.parse_args())
