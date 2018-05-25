# Libraries
import numpy as np
from datetime import datetime
from sklearn import tree

# Modules
import data
import data_reader
import process_data

class DecisionTree():
    """docstring for DecisionTree"""
    def __init__(self, filename,
                       CULTURE_SIZE_CUTOFF=0,
                       AB_CULTURE_COUNT_CUTOFF=0,
                       ESBL_AB_RESISTENCE_LIST = ["cefotaxim", "ceftazidim", "ceftriaxon"],
                       cross_validation=300,
                       testmode="cross_validation"):
        self.filename = filename
        self.data = data.Representation()
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
                                                                                                         break_amount=1,
                                                                                                         split_percentage=10)

            clf = tree.DecisionTreeClassifier()
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

            clf = tree.DecisionTreeClassifier()
            clf.fit(train_data, train_labels)

            pred_test_probabilities = clf.predict_proba(self.patient_list)

            iteration_result = []
            for row in pred_test_probabilities:
                if row[0] > row[1]:
                    iteration_result.append(["NEG", row[0]])
                else:
                    iteration_result.append(["POS", row[1]])

            predictor_result.append(iteration_result)

        # Process results
        patient_pos_probability = []
        patient_certainty = []
        for patient_index in range(len(self.patient_list)):
            pos_count = 0
            certainty_sum = 0
            for iteration_result in predictor_result:
                if iteration_result[patient_index][0] == "POS":
                    pos_count+=1
                    certainty_sum+=iteration_result[patient_index][1]
            patient_pos_probability.append(pos_count/float(len(predictor_result)))
            patient_certainty.append(certainty_sum/float(len(predictor_result)))

        return patient_pos_probability, patient_certainty

    def run(self):
        self.data.load_culture_data(self.filename)

        print("Converting data ...")
        self.data.load_filtered_esbl_patient_data()

        if self.testmode == "cross_validation":
            print("Running cross validation " + str(self.cross_validation) + " times ...")
            pos_predictor_result, neg_predictor_result = self.run_cross_validation()
           
            print("RESULTS OF TRAINING DECISION TREE:")
            print("Average accuracy of Esbl pos: " + str(np.average(pos_predictor_result)))
            print("Standard deviation of Esbl pos: " + str(np.std(pos_predictor_result)))
            print("Average accuracy of Esbl neg: " + str(np.average(neg_predictor_result)))
            print("Standard deviation of Esbl neg: " + str(np.std(neg_predictor_result)))

        elif self.testmode == "test_patient":
            self.patient_list = [self.data.esbl_pos_patient_data[3], self.data.esbl_neg_patient_data[2342]]

            temp = []
            for train_patient in self.data.esbl_pos_patient_data:
                duplicate_value = False
                for test_patient in self.patient_list:
                    if test_patient == train_patient:
                        duplicate_value = True
                if not duplicate_value: temp.append(train_patient)

            if len(temp) < len(self.data.esbl_pos_patient_data):
                print("Removed {} patient(s) from training data ...".format(len(self.data.esbl_pos_patient_data) - len(temp)))
                self.data.esbl_pos_patient_data = temp
                print("Now training on {} patients ...".format(len(self.data.esbl_pos_patient_data)))

            print("Running on test set of {} patient(s) ...".format(len(self.patient_list)))

            patient_pos_probability, patient_certainty = self.run_single_patient()

            print("RESULTS OF DECISION TREE:")
            for patient_index in range(len(patient_pos_probability)):
                print("-Patient {}".format(patient_index))
                print("Esbl pos probability: {}".format(patient_pos_probability[patient_index]))
                print("with certainty: {}".format(patient_certainty[patient_index]))