# Libraries
import numpy as np
from datetime import datetime
from sklearn.neural_network import MLPClassifier
from sklearn import tree

# Modules
import data
import data_reader
import process_data

class Perceptron():
    """docstring for Perceptron"""
    def __init__(self, filename,
                       CULTURE_SIZE_CUTOFF=0,
                       AB_CULTURE_COUNT_CUTOFF=0,
                       ESBL_AB_RESISTENCE_LIST = ["cefotaxim", "ceftazidim", "ceftriaxon"],
                       cross_validation=100,
                       testmode="cross_validation"):
        self.data = data.Representation(filename)
        self.data.set_culture_parameters(CULTURE_SIZE_CUTOFF=CULTURE_SIZE_CUTOFF,
                                         AB_CULTURE_COUNT_CUTOFF=AB_CULTURE_COUNT_CUTOFF,
                                         ESBL_AB_RESISTENCE_LIST=ESBL_AB_RESISTENCE_LIST)
        self.cross_validation = cross_validation
        self.testmode = testmode
        self.patient_list = []

    def run_cross_validation(self):
        pos_predictor_result = []
        neg_predictor_result = []
        for iteration in range(self.cross_validation):

            if iteration%10 is 0:
                print(str(iteration) + " ...")

            train_data, pos_test_data, neg_test_data, train_labels = process_data.generate_training_data(self.data,
                                                                                                         break_amount=1.2,
                                                                                                         split_percentage=20)

            clf = MLPClassifier(solver='lbfgs', alpha=1e-3, hidden_layer_sizes=(len(self.data.ab_names), 2), random_state=1)
            clf.fit(train_data, train_labels)

            pos_pred_test = clf.predict(pos_test_data)
            neg_pred_test = clf.predict(neg_test_data)

            pos_predictor_result.append(len([1 for x in pos_pred_test if x == "POS"])/float(len(pos_pred_test)))
            neg_predictor_result.append(len([1 for x in neg_pred_test if x == "NEG"])/float(len(neg_pred_test)))

        return pos_predictor_result, neg_predictor_result

    def run_single_patient(self):
        predictor_result = []

        for iteration in range(self.cross_validation):

            if iteration%10 is 0:
                print(str(iteration) + " ...")

            train_data, pos_test_data, neg_test_data, train_labels = process_data.generate_training_data(self.data,
                                                                                                         break_amount=1,
                                                                                                         split_percentage=0)

            clf = MLPClassifier(solver='lbfgs', alpha=1e-3, hidden_layer_sizes=(15, 5), random_state=1)
            clf.fit(train_data, train_labels)

            pred_test_probabilities = clf.predict_proba(self.patient_list)

            iteration_result = []
            for row in pred_test_probabilities:
                if row[0] > row[1]:
                    iteration_result.append(["NEG", row[0]])
                else:
                    iteration_result.append(["POS", row[1]])

            predictor_result.append(iteration_result)
        
        return predictor_result


    def run(self):        
        print("Converting data ...")

        self.data.load_filtered_esbl_patient_data()

        if self.testmode == "cross_validation":
            print("Running cross validation " + str(self.cross_validation) + " times ...")
            pos_predictor_result, neg_predictor_result = self.run_cross_validation()
            
            print("RESULTS OF MULTILAYER PERCEPTRON:")
            print("Average accuracy of Esbl pos: " + str(np.average(pos_predictor_result)))
            print("Standard deviation of Esbl pos: " + str(np.std(pos_predictor_result)))
            print("Average accuracy of Esbl neg: " + str(np.average(neg_predictor_result)))
            print("Standard deviation of Esbl neg: " + str(np.std(neg_predictor_result)))
        
        elif self.testmode == "test_patient":
            self.patient_list = [self.data.esbl_pos_patient_data[0]]

            temp = []
            for train_patient in self.data.esbl_pos_patient_data:
                for test_patient in self.patient_list:
                    if test_patient != train_patient:
                        temp.append(train_patient)

            if len(temp) < len(self.data.esbl_pos_patient_data):
                self.data.esbl_pos_patient_data = temp
                print("Removed {} patient(s) from training data ...".format(len(self.data.esbl_pos_patient_data) - len(temp)))

            print("Running on test set of {} patient(s) ...".format(len(self.patient_list)))
            predictor_result = self.run_single_patient()
            print("RESULTS OF MULTILAYER PERCEPTRON:")
            print(predictor_result)