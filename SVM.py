# Libraries
import numpy as np
import pandas as pd
from datetime import datetime

from sklearn import svm
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

        pos_train_data, pos_test_data = process_data.percentage_split([], [], esbl_pos_patient_data)
        if len(pos_test_data) < 1: continue
        neg_train_data, neg_test_data = process_data.percentage_split([], [], esbl_neg_patient_data)

        clf_pos = svm.OneClassSVM()
        clf_pos.fit(pos_train_data)
        clf_neg = svm.OneClassSVM()
        clf_neg.fit(neg_train_data)
        pos_pred_test = clf_pos.predict(pos_test_data)
        neg_pred_test = clf_neg.predict(pos_test_data)

        pos_predictor_result.append(len([1 for x in pos_pred_test if x>0])/float(len(pos_pred_test)))
        neg_predictor_result.append(len([1 for x in neg_pred_test if x>0])/float(len(neg_pred_test)))
    return pos_predictor_result, neg_predictor_result

def run(filename, CULTURE_SIZE_CUTOFF, AB_CULTURE_COUNT_CUTOFF, ESBL_AB_RESISTENCE_LIST, cross_validation=100):
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

    print("Running cross validation " + str(cross_validation) + " times ...")
    pos_predictor_result, neg_predictor_result = run_cross_validation(esbl_pos_patient_data, esbl_neg_patient_data, cross_validation)
    
    print("RESULTS OF TRAINING SVM:")
    print("Average accuracy of Esbl pos: " + str(np.mean(pos_predictor_result)))
    print("Standard deviation of Esbl pos: " + str(np.std(pos_predictor_result)))
    print("Average accuracy of Esbl neg: " + str(np.mean(neg_predictor_result)))
    print("Standard deviation of Esbl neg: " + str(np.std(neg_predictor_result)))