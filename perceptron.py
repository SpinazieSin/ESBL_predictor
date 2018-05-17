# Libraries
import numpy as np
import pandas as pd
from datetime import datetime

from sklearn.neural_network import MLPClassifier
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

        train_data, pos_test_data = process_data.percentage_split([], [], esbl_pos_patient_data)
        if len(pos_test_data) < 1: continue
        train_labels = ["POS" for i in range(len(train_data))]
        train_data, neg_test_data = process_data.percentage_split(train_data, [], esbl_neg_patient_data, break_limit=len(train_labels))
        for i in range(len(train_data)-len(train_labels)):
            train_labels.append("NEG")

        clf = MLPClassifier(solver='lbfgs', alpha=1e-3, hidden_layer_sizes=(32, 5), random_state=2)
        clf.fit(train_data, train_labels)
        pos_pred_test = clf.predict(pos_test_data)
        neg_pred_test = clf.predict(neg_test_data)
        pos_predictor_result.append(len([1 for x in pos_pred_test if x == "POS"])/float(len(pos_pred_test)))
        neg_predictor_result.append(len([1 for x in neg_pred_test if x == "NEG"])/float(len(neg_pred_test)))
    return pos_predictor_result, neg_predictor_result

def run(filename, CULTURE_SIZE_CUTOFF, AB_CULTURE_COUNT_CUTOFF, ESBL_AB_RESISTENCE_LIST, cross_validation=50):
    csv_data = data_reader.read_csv(filename)
    
    print("Converting data ...")
    id_dict = csv_data[0]
    ab_dict = csv_data[1]
    ab_names = ab_dict.keys()

    esbl_pos_patient_data, esbl_neg_patient_data = process_data.generate_esbl_patient_data(id_dict, ab_dict, CULTURE_SIZE_CUTOFF, ESBL_AB_RESISTENCE_LIST)

    esbl_pos_ab_count = process_data.count_ab(esbl_pos_patient_data, ab_names)
    esbl_neg_ab_count = process_data.count_ab(esbl_neg_patient_data, ab_names)
    relevant_ab_list = process_data.relevant_ab(esbl_pos_ab_count, ab_dict, ab_names, AB_CULTURE_COUNT_CUTOFF)

    esbl_pos_patient_data = process_data.filter_ab(esbl_pos_patient_data, relevant_ab_list, ab_dict, None, numeric=True)
    esbl_neg_patient_data = process_data.filter_ab(esbl_neg_patient_data, relevant_ab_list, ab_dict, None, numeric=True)

    print("Running cross validation " + str(cross_validation) + " times ...")
    pos_predictor_result, neg_predictor_result = run_cross_validation(esbl_pos_patient_data, esbl_neg_patient_data, cross_validation)

    print("RESULTS OF TRAINING DECISION TREE:")
    print("Average accuracy of Esbl pos: " + str(np.mean(pos_predictor_result)))
    print("Standard deviation of Esbl pos: " + str(np.std(pos_predictor_result)))
    print("Average accuracy of Esbl neg: " + str(np.mean(neg_predictor_result)))
    print("Standard deviation of Esbl neg: " + str(np.std(neg_predictor_result)))