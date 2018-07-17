# Libraries
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from sklearn.neural_network import MLPClassifier
from random import randint
from time import time

# Modules
import data
import data_reader
import process_data

class Perceptron():
    def __init__(self, filename,
                       CULTURE_SIZE_CUTOFF=0,
                       AB_CULTURE_COUNT_CUTOFF=0,
                       ESBL_AB_RESISTANCE_LIST = ["cefotaxim", "ceftazidim", "ceftriaxon"],
                       RELEVANT_MO_LIST=[],
                       testmode="cross_validation",
                       analysis_type="culture",
                       medication_file=None,
                       cross_validation=3000):
        self.filename = filename
        self.data = data.Representation()
        self.data.set_culture_parameters(CULTURE_SIZE_CUTOFF=CULTURE_SIZE_CUTOFF,
                                         AB_CULTURE_COUNT_CUTOFF=AB_CULTURE_COUNT_CUTOFF,
                                         ESBL_AB_RESISTANCE_LIST=ESBL_AB_RESISTANCE_LIST,
                                         RELEVANT_MO_LIST=RELEVANT_MO_LIST,)
        self.testmode = testmode
        self.analysis_type = analysis_type
        self.medication_file = medication_file

        self.cross_validation = cross_validation
        self.layers = [8, 2]

        self.iteration_timer_estimate = 0

    def run_cross_validation(self):
        pos_predictor_result = []
        neg_predictor_result = []

        self.iteration_timer_estimate = 0

        importance = np.zeros(len(self.data.relevant_ab_names))
        
        for iteration in range(self.cross_validation):

            # Start of timer
            iteration_start = time()

            if iteration%50 is 0:
                print(str(iteration) + " ...")
                print("Time left: {}".format((self.iteration_timer_estimate/(iteration+1))*(self.cross_validation-iteration)))

            # Split dataset into train and test
            train_data, pos_test_data, neg_test_data, train_labels = process_data.generate_training_data(self.pos_training_data,
                                                                                                         self.neg_training_data,
                                                                                                         break_amount=1,
                                                                                                         split_percentage=10)

            # Train percpetron on training set
            clf = MLPClassifier(solver='lbfgs', alpha=1e-3, hidden_layer_sizes=(self.layers[0], self.layers[1]), random_state=0)
            clf.fit(train_data, train_labels)

            # Test the classification on a test and control set
            pos_pred_test = clf.predict(pos_test_data)
            neg_pred_test = clf.predict(neg_test_data)

            # Calculate accuracy and store it
            pos_predictor_result.append(len([1 for x in pos_pred_test if x == "POS"])/float(len(pos_pred_test)))
            neg_predictor_result.append(len([1 for x in neg_pred_test if x == "NEG"])/float(len(neg_pred_test)))

            self.iteration_timer_estimate+=(time()-iteration_start)

        return pos_predictor_result, neg_predictor_result

    def run_single_patient(self):
        predictor_result = []

        for iteration in range(self.cross_validation):

            if iteration%20 is 0:
                print(str(iteration) + " ...")

            # Create training set (split percentage is 0)
            train_data, pos_test_data, neg_test_data, train_labels = process_data.generate_training_data(self.pos_training_data,
                                                                                                         self.neg_training_data,
                                                                                                         break_amount=1,
                                                                                                         split_percentage=0)

            # Train a perceptron classifier
            clf = MLPClassifier(solver='lbfgs', alpha=1e-3, hidden_layer_sizes=(self.layers[0], self.layers[1]), random_state=0)
            clf.fit(train_data, train_labels)

            # Output prediction probabilities
            pred_test_probabilities = clf.predict_proba(self.patient_list)

            # Append labels
            iteration_result = []
            for row in pred_test_probabilities:
                if row[0] > row[1]:
                    iteration_result.append(["NEG", row[0]])
                else:
                    iteration_result.append(["POS", row[1]])

            predictor_result.append(iteration_result)

        # Process results, sums over all suscebtible and certainty values per patient
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
        
        if self.testmode == "date" and self.analysis_type == "medication":
            if not self.medication_file:
                print("Medication file not set")
                return
            self.data.load_patient_data(self.filename, self.medication_file)
 
            max_date = 15
            min_date = 1
            results = [[0 for _ in range(max_date)] for _ in range(min_date)]
            score = []
            sensitivity = []
            specificity = []
            sensitivity_variance = []
            specificity_variance = []
            x_axis = []

            print("Beginning date iteration")
            for min_month in range(0, min_date):
                
                for max_month in range(0, max_date):
                    if max_month < min_month: continue

                    self.data.date_range = [ 1 + min_month*30 , 2 + max_month*30 ]
                    print("New range - {}".format(self.data.date_range))
                    self.data.generate_patient_data()
                    self.pos_training_data, self.neg_training_data = process_data.vectorize_medication_data(self.data.patient_dict)

                    pos_predictor_result, neg_predictor_result = self.run_cross_validation()

                    print("Average accuracy of Esbl pos: " + str(np.average(pos_predictor_result)))
                    print("Average accuracy of Esbl neg: " + str(np.average(neg_predictor_result)))

                    # results[min_month][max_month] = (np.average(pos_predictor_result)*np.average(neg_predictor_result)*len(self.pos_training_data))/(np.std(pos_predictor_result)*np.std(neg_predictor_result)+1)
                    results[min_month][max_month] = (np.average(pos_predictor_result)*np.average(neg_predictor_result))/(np.std(pos_predictor_result)*np.std(neg_predictor_result)+1)
                    
                    # Graph results
                    score.append((np.average(pos_predictor_result)*np.average(neg_predictor_result))/(np.std(pos_predictor_result)*np.std(neg_predictor_result)+1))
                    sensitivity.append(np.mean(pos_predictor_result))
                    specificity.append(np.mean(neg_predictor_result))
                    sensitivity_variance.append(np.std(pos_predictor_result))
                    specificity_variance.append(np.std(neg_predictor_result))
                    x_axis.append(max_month)

            use_graph = True
            use_grid = False
            if use_graph:
                print(score)
                def sliding_mean(data_array, window=5):  
                    data_array = array(data_array)  
                    new_list = []  
                    for i in range(len(data_array)):  
                        indices = range(max(i - window + 1, 0),  
                                        min(i + window + 1, len(data_array)))  
                        avg = 0  
                        for j in indices:  
                            avg += data_array[j]  
                        avg /= float(len(indices))  
                        new_list.append(avg)  
                          
                    return array(new_list)  
 
                plt.figure(figsize=(12, 9))  
                  
                ax = plt.subplot(111)  
                ax.spines["top"].set_visible(False)  
                ax.spines["right"].set_visible(False)  
                  
                ax.get_xaxis().tick_bottom()  
                ax.get_yaxis().tick_left()  

                plt.ylim(0, 1)  

                plt.xticks(fontsize=14)  
                plt.yticks(fontsize=14)  

                plt.ylabel("Score", fontsize=16)  
                  

                plt.plot(x_axis, sensitivity, color="blue", lw=2)
                plt.plot(x_axis, specificity, color="green", lw=2)
                plt.plot(x_axis, score, color="red", lw=2)
                plt.title("", fontsize=22)  
                plt.xlabel("\nAge of medication", fontsize=10)  
                  
                plt.show() 

            if use_grid:
                best_score = np.max(results)
                for layer_1 in range(len(results)):
                    for layer_2 in range(len(results[layer_1])):
                        if results[layer_1][layer_2] == best_score:
                            print("\n")
                            print("Best date: {},{}".format((layer_1+1)*2, layer_2))
                            print("With Score: {}".format(results[layer_1][layer_2]))

                results/=np.max(best_score)
                fig = plt.figure(figsize=(8, 4))

                ax = fig.add_subplot(111)

                ax.set_title('colorMap')
                plt.imshow(results)
                ax.set_aspect('equal')

                cax = fig.add_axes([0.12, 0.1, 0.78, 0.8])
                cax.get_xaxis().set_visible(False)
                cax.get_yaxis().set_visible(False)
                cax.patch.set_alpha(0)
                cax.set_frame_on(False)
                plt.show()
            return

        elif self.testmode == "date":
            self.data.load_culture_data(self.filename)

            max_date = 14
            min_date = 1
            results = [[0 for _ in range(max_date)] for _ in range(min_date)]
            score = []
            sensitivity = []
            specificity = []
            sensitivity_variance = []
            specificity_variance = []
            x_axis = []

            print("Beginning date iteration")
            for min_month in range(0, min_date):
                
                for max_month in range(0, max_date):
                    if max_month < min_month: continue

                    self.data.date_range = [ 2 + min_month*30 , 3 + max_month*30 ]
                    print("New range - {}".format(self.data.date_range))

                    self.data.load_filtered_esbl_patient_data()

                    self.pos_training_data = self.data.esbl_pos_patient_data
                    if len(self.pos_training_data) < 10: continue
                    self.neg_training_data = self.data.esbl_neg_patient_data


                    pos_predictor_result, neg_predictor_result = self.run_cross_validation()

                    print("Average accuracy of Esbl pos: " + str(np.average(pos_predictor_result)))
                    print("Average accuracy of Esbl neg: " + str(np.average(neg_predictor_result)))

                    # Grid results
                    # results[min_month][max_month] = (np.average(pos_predictor_result)*np.average(neg_predictor_result)*len(self.pos_training_data))/(np.std(pos_predictor_result)*np.std(neg_predictor_result)+1)
                    results[min_month][max_month] = (np.average(pos_predictor_result)*np.average(neg_predictor_result))/(np.std(pos_predictor_result)*np.std(neg_predictor_result)+1)

                    # Graph results
                    score.append((np.average(pos_predictor_result)*np.average(neg_predictor_result))/(np.std(pos_predictor_result)*np.std(neg_predictor_result)+1))
                    sensitivity.append(np.mean(pos_predictor_result))
                    specificity.append(np.mean(neg_predictor_result))
                    sensitivity_variance.append(np.std(pos_predictor_result))
                    specificity_variance.append(np.std(neg_predictor_result))
                    x_axis.append(max_month)

            use_graph = True
            use_grid = False
            if use_graph:
                print(score)
                def sliding_mean(data_array, window=5):  
                    data_array = array(data_array)  
                    new_list = []  
                    for i in range(len(data_array)):  
                        indices = range(max(i - window + 1, 0),  
                                        min(i + window + 1, len(data_array)))  
                        avg = 0  
                        for j in indices:  
                            avg += data_array[j]  
                        avg /= float(len(indices))  
                        new_list.append(avg)  
                          
                    return array(new_list)  
 
                plt.figure(figsize=(12, 9))  
                  
                ax = plt.subplot(111)  
                ax.spines["top"].set_visible(False)  
                ax.spines["right"].set_visible(False)  
                  
                ax.get_xaxis().tick_bottom()  
                ax.get_yaxis().tick_left()  

                plt.ylim(0, 1)  

                plt.xticks(fontsize=14)  
                plt.yticks(fontsize=14)  

                plt.ylabel("Score", fontsize=16)  
                  

                plt.plot(x_axis, sensitivity, color="blue", lw=2)
                plt.plot(x_axis, specificity, color="green", lw=2)
                plt.plot(x_axis, score, color="red", lw=2)
                plt.title("", fontsize=22)  
                plt.xlabel("\nAge of cultures", fontsize=10)  
                  
                plt.show()

            if use_grid:
                best_score = np.max(results)
                for layer_1 in range(len(results)):
                    for layer_2 in range(len(results[layer_1])):
                        if results[layer_1][layer_2] == best_score:
                            print("\n")
                            print("Best date: {},{}".format((layer_1+1)*2, layer_2))
                            print("With Score: {}".format(results[layer_1][layer_2]))

                results/=np.max(best_score)
                fig = plt.figure(figsize=(8, 4))

                ax = fig.add_subplot(111)

                ax.set_title('colorMap')
                plt.imshow(results)
                ax.set_aspect('equal')

                cax = fig.add_axes([0.12, 0.1, 0.78, 0.8])
                cax.get_xaxis().set_visible(False)
                cax.get_yaxis().set_visible(False)
                cax.patch.set_alpha(0)
                cax.set_frame_on(False)
                plt.show()
            return

        if self.analysis_type=="culture" and self.testmode != "date":
            self.data.load_culture_data(self.filename)

            print("Converting data ...")
            self.data.load_filtered_esbl_patient_data()
            self.pos_training_data = self.data.esbl_pos_patient_data
            self.neg_training_data = self.data.esbl_neg_patient_data

        elif self.analysis_type=="medication" and self.testmode != "date":
            if not self.medication_file:
                print("Medication file not set")
                return
            self.data.load_patient_data(self.filename, self.medication_file)
            self.pos_training_data, self.neg_training_data = process_data.vectorize_medication_data(self.data.patient_dict)
        else:
            print("Invalid analysis type")
            return

        if self.testmode == "cross_validation":
            print("Running cross validation " + str(self.cross_validation) + " times ...")

            pos_predictor_result, neg_predictor_result = self.run_cross_validation()
            
            print("RESULTS OF MULTILAYER PERCEPTRON:")
            print("Average accuracy of Esbl pos: " + str(np.average(pos_predictor_result)))
            print("Standard deviation of Esbl pos: " + str(np.std(pos_predictor_result)))
            print("Average accuracy of Esbl neg: " + str(np.average(neg_predictor_result)))
            print("Standard deviation of Esbl neg: " + str(np.std(neg_predictor_result)))
        
        elif self.testmode == "test_patient":

            for patient in self.data.esbl_pos_patient_data:
                print("____________________________")
                self.patient_list = [patient, self.data.esbl_neg_patient_data[randint(0, len(self.data.esbl_neg_patient_data)-1)]]
                
                self.pos_training_data = self.data.esbl_pos_patient_data
                
                temp = []
                for train_patient in self.pos_training_data:
                    duplicate_value = False
                    for test_patient in self.patient_list:
                        if test_patient == train_patient:
                            duplicate_value = True
                    if not duplicate_value: temp.append(train_patient)

                if len(temp) < len(self.pos_training_data):
                    print("Removed {} patient(s) from training data ...".format(len(self.pos_training_data) - len(temp)))
                    self.pos_training_data = temp
                    print("Now training on {} patients ...".format(len(self.pos_training_data)))

                print("Running on test set of {} patient(s) ...".format(len(self.patient_list)))

                patient_pos_probability, patient_certainty = self.run_single_patient()

                print("RESULTS OF MULTILAYER PERCEPTRON:")
                for patient_index in range(len(patient_pos_probability)):
                    print("-Patient {}".format(patient_index))
                    print("Was resistent to: {}".format([self.data.relevant_ab_names[ab_index] for ab_index in range(len(self.patient_list[patient_index])) if self.patient_list[patient_index][ab_index] < 0] ))
                    print("Pos probability: {}".format(patient_pos_probability[patient_index]))
                    print("with certainty: {}".format(patient_certainty[patient_index]))


        elif self.testmode == "perceptron_optimization":
            max_layer = 15
            results = [[0 for _ in range(max_layer-1)] for _ in range(max_layer-1)]

            print("Beginning date iteration")
            for layer_1 in range(1, max_layer):
                
                for layer_2 in range(1, max_layer):

                    self.layers = [ 2*layer_1, 2*layer_2 ]

                    print("New layers - {}".format(self.layers))

                    pos_predictor_result, neg_predictor_result = self.run_cross_validation()
                    print(pos_predictor_result)
                    print(neg_predictor_result)
                    print("Average accuracy of pos: " + str(np.average(pos_predictor_result)))
                    print("Average accuracy of neg: " + str(np.average(neg_predictor_result)))

                    results[layer_1-1][layer_2-1] = (np.average(pos_predictor_result)*np.average(neg_predictor_result))/(np.std(pos_predictor_result)*np.std(neg_predictor_result)+1)
            
            for layer_1 in range(len(results)):
                for layer_2 in range(len(results[layer_1])):
                    if results[layer_1][layer_2] == np.max(results):
                        print("Best layers: {},{}".format((layer_1+1)*2, layer_2))
            results/=np.max(results)

            fig = plt.figure(figsize=(6, 3.2))

            ax = fig.add_subplot(111)

            ax.set_title('colorMap')
            plt.imshow(results)
            ax.set_aspect('equal')

            cax = fig.add_axes([0.12, 0.1, 0.78, 0.8])
            cax.get_xaxis().set_visible(False)
            cax.get_yaxis().set_visible(False)
            cax.patch.set_alpha(0)
            cax.set_frame_on(False)
            plt.colorbar(orientation='vertical')
            plt.show()
            return