# Libraries
import numpy as np
import pandas as pd
from datetime import datetime

from sklearn import tree
import matplotlib.pyplot as plt
import matplotlib.font_manager
# Modules
import data_reader
import process_data
import data_reader


def run_cross_validation(esbl_pos_patient_data, esbl_neg_patient_data, cross_validation):
    pos_predictor_result = []
    neg_predictor_result = []
    for iteration in range(cross_validation):

        if iteration%10 is 0:
            print(str(iteration) + " ...")

        train_data, pos_test_data, neg_test_data, train_labels = process_data.generate_training_data(esbl_pos_patient_data,
                                                                                                     esbl_neg_patient_data,
                                                                                                     break_amount=2)
        clf = tree.DecisionTreeClassifier()
        clf.fit(train_data, train_labels)

        pos_pred_test = clf.predict(pos_test_data)
        neg_pred_test = clf.predict(neg_test_data)

        pos_predictor_result.append(len([1 for x in pos_pred_test if x == "POS"])/float(len(pos_pred_test)))
        neg_predictor_result.append(len([1 for x in neg_pred_test if x == "NEG"])/float(len(neg_pred_test)))

    return pos_predictor_result, neg_predictor_result

def run_single_patient(esbl_pos_patient_data, esbl_neg_patient_data, iterations, patient_list):
    predictor_result = []

    for iteration in range(iterations):

        if iteration%10 is 0:
            print(str(iteration) + " ...")

        train_data, pos_test_data, neg_test_data, train_labels = process_data.generate_training_data(esbl_pos_patient_data,
                                                                                                     esbl_neg_patient_data,
                                                                                                     break_amount=1,
                                                                                                     split_percentage=0)


        clf = tree.DecisionTreeClassifier()
        clf.fit(train_data, train_labels)

        pred_test_probabilities = clf.predict_proba(patient_list)

        iteration_result = []
        for row in pred_test_probabilities:
            if row[0] > row[1]:
                iteration_result.append(["NEG", row[0]])
            else:
                iteration_result.append(["POS", row[1]])

        predictor_result.append(iteration_result)
    
    return predictor_result

def run(filename, CULTURE_SIZE_CUTOFF, AB_CULTURE_COUNT_CUTOFF, ESBL_AB_RESISTENCE_LIST, cross_validation=100, testmode="cross_validation"):
    csv_data = data_reader.read_csv(filename)
    
    print("Converting data ...")
    id_dict = csv_data[0]
    ab_dict = csv_data[1]
    ab_names = ab_dict.keys()

    esbl_pos_patient_data, esbl_neg_patient_data = process_data.generate_data(patient_data=id_dict,
                                                                                           ab_data=ab_dict,
                                                                                           ab_names=ab_names,
                                                                                           AB_CULTURE_COUNT_CUTOFF=AB_CULTURE_COUNT_CUTOFF,
                                                                                           CULTURE_SIZE_CUTOFF=CULTURE_SIZE_CUTOFF,
                                                                                           ESBL_AB_RESISTENCE_LIST=ESBL_AB_RESISTENCE_LIST,
                                                                                           esbl_result_format=None,
                                                                                           numeric=True)

    if testmode == "cross_validation":
        print("Running cross validation " + str(cross_validation) + " times ...")
        pos_predictor_result, neg_predictor_result = run_cross_validation(esbl_pos_patient_data, esbl_neg_patient_data, cross_validation)

        print("RESULTS OF TRAINING DECISION TREE:")
        print("Average accuracy of Esbl pos: " + str(np.mean(pos_predictor_result)))
        print("Standard deviation of Esbl pos: " + str(np.std(pos_predictor_result)))
        print("Average accuracy of Esbl neg: " + str(np.mean(neg_predictor_result)))
        print("Standard deviation of Esbl neg: " + str(np.std(neg_predictor_result)))

    elif testmode == "test_patient":
        patient_list = [esbl_neg_patient_data[3767]]
        temp = []
        for train_patient in esbl_pos_patient_data:
            for test_patient in patient_list:
                if test_patient != train_patient:
                    temp.append(train_patient)

        if len(temp) < len(esbl_pos_patient_data):
            print("Removed {} patient(s) from training data ...".format(len(esbl_pos_patient_data) - len(temp)))

        print("Running on test set of {} patient(s) ...".format(len(patient_list)))
        predictor_result = run_single_patient(temp, esbl_neg_patient_data, cross_validation, patient_list)
        print("RESULTS OF TRAINING DECISION TREE:")
        print(predictor_result)