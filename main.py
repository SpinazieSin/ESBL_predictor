# Libraries
import argparse

# Modules
import average
import SVM
import decision_tree
import perceptron

# Defined constants
ESBL_AB_RESISTANCE_LIST = ["cefotaxim", "ceftazidim"] # 3rde generatie cefalosporine resistentie
CULTURE_SIZE_CUTOFF = 5
AB_CULTURE_COUNT_CUTOFF = 20

def main(args):
    filename = args.file
    method = args.method

    if args.average:
        classifier = average.Average(filename, args.culture_cutoff, args.ab_count_cutoff, ESBL_AB_RESISTANCE_LIST)
    # elif args.svm:
    #     SVM.run(filename, args.culture_cutoff, args.ab_count_cutoff, ESBL_AB_RESISTANCE_LIST)
    elif args.tree:
        classifier = decision_tree.DecisionTree(filename, args.culture_cutoff, args.ab_count_cutoff, ESBL_AB_RESISTANCE_LIST, testmode=method)
    elif args.perceptron:
        classifier = perceptron.Perceptron(filename, args.culture_cutoff, args.ab_count_cutoff, ESBL_AB_RESISTANCE_LIST, testmode=method)
    
    classifier.run()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-file', type=str, help='File to use for analysis', default="data/sample.csv")
    parser.add_argument('-average', action='store_true', help='Train and print results from baseline analysis method')
    # parser.add_argument('-svm', action='store_true', help='Train and print results from SVM method')
    parser.add_argument('-tree', action='store_true', help='Train and print results from decision tree method')
    parser.add_argument('-perceptron', action='store_true', help='Train and print results from multilayer perceptron neural network')
    parser.add_argument('-culture_cutoff', type=int, help='Determines the minimum amount of cultures per patient', default=5)
    parser.add_argument('-ab_count_cutoff', type=int, help='In what percentage of cultures should an AB occur in order to be relevant?', default=20)
    parser.add_argument('-method', type=str, help='Either test_patient or cross_validation', default='cross_validation')
    main(parser.parse_args())
