# Libraries
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from sklearn import tree
from random import randint

# Modules
import data
import data_reader
import process_data

class DecisionTree():
    def __init__(self, filename,
                       CULTURE_SIZE_CUTOFF=0,
                       AB_CULTURE_COUNT_CUTOFF=0,
                       ESBL_AB_RESISTANCE_LIST = ["cefotaxim", "ceftazidim", "ceftriaxon"],
                       testmode="cross_validation",
                       analysis_type="culture",
                       medication_file=None,
                       cross_validation=2000):
        self.filename = filename
        self.data = data.Representation()
        self.data.set_culture_parameters(CULTURE_SIZE_CUTOFF=CULTURE_SIZE_CUTOFF,
                                         AB_CULTURE_COUNT_CUTOFF=AB_CULTURE_COUNT_CUTOFF,
                                         ESBL_AB_RESISTANCE_LIST=ESBL_AB_RESISTANCE_LIST)
        self.testmode = testmode
        self.analysis_type=analysis_type
        self.cross_validation = cross_validation

    def run_cross_validation(self):
        pos_predictor_result = []
        neg_predictor_result = []

        importance = np.zeros(len(self.data.relevant_ab_names))
        for iteration in range(self.cross_validation):

            if iteration%500 is 0:
                print(str(iteration) + " ...")

            train_data, pos_test_data, neg_test_data, train_labels = process_data.generate_training_data(self.pos_training_data,
                                                                                                         self.neg_training_data,
                                                                                                         break_amount=1,
                                                                                                         split_percentage=10)

            clf = tree.DecisionTreeClassifier()
            clf.fit(train_data, train_labels)
            importance+=clf.feature_importances_

            pos_pred_test = clf.predict(pos_test_data)
            neg_pred_test = clf.predict(neg_test_data)

            pos_predictor_result.append(len([1 for x in pos_pred_test if x == "POS"])/float(len(pos_pred_test)))
            neg_predictor_result.append(len([1 for x in neg_pred_test if x == "NEG"])/float(len(neg_pred_test)))

        print("")
        results = []
        for i in range(len(self.data.relevant_ab_names)):
            results.append([importance[i]/self.cross_validation, self.data.relevant_ab_names[i]])
        results.sort(reverse=True)
        for i in range(len(self.data.relevant_ab_names)):
            print("{}: {}".format(results[i][1], results[i][0]))
        print("")
        return pos_predictor_result, neg_predictor_result

    def run_single_patient(self):
        predictor_result = []

        for iteration in range(self.cross_validation):

            if iteration%500 is 0:
                print(str(iteration) + " ...")

            train_data, pos_test_data, neg_test_data, train_labels = process_data.generate_training_data(self.pos_training_data,
                                                                                                         self.neg_training_data,
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

        if self.testmode == "date":
            self.data.load_culture_data(self.filename)

            max_date = 5
            pos_results = [[0 for _ in range(max_date)] for _ in range(max_date)]
            neg_results = [[0 for _ in range(max_date)] for _ in range(max_date)]
            print("Beginning date iteration")
            for min_month in range(0, max_date):
                
                for max_month in range(0, max_date):
                    if max_month < min_month: continue

                    date_range=[ 1 + min_month*30 , 1 + max_month*30 ]
                    print("New range - {}".format(date_range))

                    self.data.load_filtered_esbl_patient_data(date_range=date_range)
                    self.pos_training_data = self.data.esbl_pos_patient_data
                    if len(self.pos_training_data) < 2: continue
                    self.neg_training_data = self.data.esbl_neg_patient_data
                    
                    pos_predictor_result, neg_predictor_result = self.run_cross_validation()
                    
                    print("Average accuracy of Esbl pos: " + str(np.average(pos_predictor_result)))
                    print("Average accuracy of Esbl neg: " + str(np.average(neg_predictor_result)))

                    pos_results[min_month][max_month] = np.average(pos_predictor_result)
                    neg_results[min_month][max_month] = np.average(neg_predictor_result)
            
            print(pos_results)
            fig = plt.figure(figsize=(6, 3.2))

            ax = fig.add_subplot(111)

            ax.set_title('colorMap')
            plt.imshow(pos_results)
            ax.set_aspect('equal')

            cax = fig.add_axes([0.12, 0.1, 0.78, 0.8])
            cax.get_xaxis().set_visible(False)
            cax.get_yaxis().set_visible(False)
            cax.patch.set_alpha(0)
            cax.set_frame_on(False)
            plt.colorbar(orientation='vertical')
            plt.show()
            return

        if self.analysis_type=="culture":
            self.data.load_culture_data(self.filename)
            
            print("Converting data ...")
            self.data.load_filtered_esbl_patient_data()
            self.pos_training_data = self.data.esbl_pos_patient_data
            self.neg_training_data = self.data.esbl_neg_patient_data

        else:
            if not medication_file:
                print("Medication file not set")
                return
            self.data.load_patient_data("data/ecoli_resistentie_bepaling.csv", "data/medicatie_toediening.csv")
            self.pos_training_data, self.neg_training_data = process_data.vectorize_medication_data(self.data.patient_dict)

        if self.testmode == "cross_validation":
            print("Running cross validation " + str(self.cross_validation) + " times ...")

            pos_predictor_result, neg_predictor_result = self.run_cross_validation()
           
            print("RESULTS OF TRAINING DECISION TREE:")
            print("Average accuracy of Esbl pos: " + str(np.mean(pos_predictor_result)))
            print("Standard deviation of Esbl pos: " + str(np.std(pos_predictor_result)))
            print("Average accuracy of Esbl neg: " + str(np.mean(neg_predictor_result)))
            print("Standard deviation of Esbl neg: " + str(np.std(neg_predictor_result)))

        elif self.testmode == "test_patient":

            for patient in self.data.esbl_pos_patient_data:
                print("____________________________")
                self.patient_list = [patient, self.data.esbl_neg_patient_data[randint(0, len(self.data.esbl_neg_patient_data)-1)]]

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
                    print("Was resistent to: {}".format([self.data.relevant_ab_names[ab_index] for ab_index in range(len(self.patient_list[patient_index])) if self.patient_list[patient_index][ab_index] < 0] ))
                    print("Pos probability: {}".format(patient_pos_probability[patient_index]))
                    print("with certainty: {}".format(patient_certainty[patient_index]))