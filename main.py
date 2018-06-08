# Libraries
import argparse

# Modules
import average
import SVM
import decision_tree
import perceptron

# Defined constants
ESBL_AB_RESISTANCE_LIST = ["cefotaxim", "ceftazidim", "ceftriaxon"] # 3rde generatie cefalosporine resistentie
# RELEVANT_MO_LIST = ["Escherichia coli"]
RELEVANT_MO_LIST = []
CULTURE_SIZE_CUTOFF = 10
AB_CULTURE_COUNT_CUTOFF = 50

def main(args):
    filename = args.file
    method = args.method

    print("-STARTING-\n")
    print("Using")
    print("Ab cutoff: {}".format(args.ab_count_cutoff))
    print("Culture cutoff: {}".format(args.culture_cutoff))
    print("Method: {}".format(method))
    print("Analysis type: {}".format(args.analysis_type))
    print()

    if args.average:
        classifier = average.Average(filename,
                                     CULTURE_SIZE_CUTOFF=args.culture_cutoff,
                                     AB_CULTURE_COUNT_CUTOFF=args.ab_count_cutoff,
                                     ESBL_AB_RESISTENCE_LIST=ESBL_AB_RESISTANCE_LIST,
                                     RELEVANT_MO_LIST=RELEVANT_MO_LIST)
    # elif args.svm:
    #     SVM.run(filename, args.culture_cutoff, args.ab_count_cutoff, ESBL_AB_RESISTANCE_LIST)
    elif args.tree:
        classifier = decision_tree.DecisionTree(filename,
                                                CULTURE_SIZE_CUTOFF=args.culture_cutoff,
                                                AB_CULTURE_COUNT_CUTOFF=args.ab_count_cutoff,
                                                ESBL_AB_RESISTANCE_LIST=ESBL_AB_RESISTANCE_LIST,
                                                RELEVANT_MO_LIST=RELEVANT_MO_LIST,
                                                testmode=method,
                                                analysis_type=args.analysis_type,
                                                medication_file=args.medication_file)
    elif args.perceptron:
        classifier = perceptron.Perceptron(filename,
                                           CULTURE_SIZE_CUTOFF=args.culture_cutoff,
                                           AB_CULTURE_COUNT_CUTOFF=args.ab_count_cutoff,
                                           ESBL_AB_RESISTANCE_LIST=ESBL_AB_RESISTANCE_LIST,
                                           RELEVANT_MO_LIST=RELEVANT_MO_LIST,
                                           testmode=method,
                                           analysis_type=args.analysis_type,
                                           medication_file=args.medication_file)
    
    classifier.run()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-file', type=str, help='File to use for culture analysis', default="data/culture_sample.csv")
    parser.add_argument('-average', action='store_true', help='Train and print results from baseline analysis method')
    # parser.add_argument('-svm', action='store_true', help='Train and print results from SVM method')
    parser.add_argument('-tree', action='store_true', help='Train and print results from decision tree method')
    parser.add_argument('-perceptron', action='store_true', help='Train and print results from multilayer perceptron neural network')
    parser.add_argument('-culture_cutoff', type=int, help='Determines the minimum amount of cultures per patient', default=10)
    parser.add_argument('-ab_count_cutoff', type=int, help='In what percentage of cultures should an AB occur in order to be relevant?', default=50)
    parser.add_argument('-method', type=str, help='Either test_patient, cross_validation or date', default='cross_validation')
    parser.add_argument('-analysis_type', type=str, help='Either culture or medication', default='culture')
    parser.add_argument('-medication_file', type=str, help='If analysis_type is medication, this file needs to be set', default="data/medication_sample.csv")
    main(parser.parse_args())
