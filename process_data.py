# Libraries
import numpy as np
import itertools as it
from random import randint
from math import floor

# Modules
import data

def count_ab(esbl_patient_data, ab_names):
    ab_count = np.zeros((len(ab_names), 2))
    for ab_vector in esbl_patient_data:
        if np.count_nonzero(ab_vector) > 0:
            for ab_index in range(len(ab_vector)):
                if ab_vector[ab_index] is 'S':
                    # Suscebtible count
                    ab_count[ab_index][0]+=1
                    # Resistent count
                    ab_count[ab_index][1]+=1
                elif ab_vector[ab_index] is 'R':
                    # Reistent count
                    ab_count[ab_index][1]+=1
    return ab_count

def relevant_ab(esbl_pos_ab_count, ab_dict, ab_names, RELATIVE_AB_CULTURE_COUNT_CUTOFF):
    relevant_ab_list = []
    for ab in ab_names:
        ab_index = ab_dict[ab]
        if esbl_pos_ab_count[ab_index][1] > RELATIVE_AB_CULTURE_COUNT_CUTOFF:
            relevant_ab_list.append(ab)
    return relevant_ab_list

def filter_ab(patient_data, relevant_ab_list, ab_dict, esbl_result, numeric):
    esbl_patient_data = []
    for row in patient_data:
        patient_ab = []
        for ab in relevant_ab_list:
            ab_result = row[ab_dict[ab]]
            if numeric:
                value = 0
                if ab_result is 'S':
                    value = 1
                elif ab_result is 'R':
                    value = -1
                patient_ab.append(value)
            else:
                patient_ab.append(row[ab_dict[ab]])
        if esbl_result is not None:
            patient_ab.append(esbl_result)
        esbl_patient_data.append(patient_ab)
    return esbl_patient_data

def percentage_split(train_data, test_data, esbl_patient_data, split_percentage, sample_count=-1, random_sampling=False):
    break_count = 0
    if sample_count > 0 and random_sampling:
        data_length = len(esbl_patient_data)-1
        for i in range(sample_count):
            if randint(0, 100) < split_percentage:
                test_data.append(esbl_patient_data[randint(0, data_length)])
            else:
                train_data.append(esbl_patient_data[randint(0, data_length)])
    else:
        for data_row in esbl_patient_data:
            if break_count == sample_count: break
            break_count+=1
            if randint(0, 100) < split_percentage:
                test_data.append(data_row)
            else:
                train_data.append(data_row)
    return train_data, test_data 

def find_esbl_pos_day(patient_data, ESBL_AB_RESISTENCE_LIST):
    for culture in patient_data:
        ab = culture[1]
        result = culture[2]
        bepaling = culture[3]
        if (ab in ESBL_AB_RESISTENCE_LIST and result is not 'S' and "M_banaal_BK" in bepaling):
            return culture[0]
    return None


def generate_esbl_patient_data(id_dict, ab_dict, CULTURE_SIZE_CUTOFF, ESBL_AB_RESISTENCE_LIST):
    esbl_pos_patient_data = []
    esbl_neg_patient_data = []
    ab_length = len(ab_dict.keys())
    
    for patient in id_dict.keys():

        esbl_found = False
        patient_data = id_dict[patient]
        if len(patient_data) < CULTURE_SIZE_CUTOFF:
            continue
        
        patient_data.sort()
        ab_vector = [None for x in range(ab_length)]

        esbl_pos_day = find_esbl_pos_day(patient_data, ESBL_AB_RESISTENCE_LIST)
        previous_culture_day = None

        for culture in patient_data:

            culture_day = culture[0]
            ab = culture[1]
            result = culture[2]
            if culture_day == esbl_pos_day:
                esbl_found = True
                # Remove this break if you want to include the cultures of the day
                # where the patient had the ESBL producing microorganism.
                break

            if (previous_culture_day != culture_day):
                previous_culture_day = culture_day
            skip = False
            if (previous_culture_day is not None and esbl_pos_day is not None):
                total_days= abs((previous_culture_day-esbl_pos_day).days)
                if (total_days < 2 or total_days > 360):
                    pass
                    skip = True

            if not skip:
                if result is 'S':
                    ab_vector[ab_dict[ab]] = 'S'
                else:
                    ab_vector[ab_dict[ab]] = 'R'

            
        if np.count_nonzero(ab_vector) > 0:
            if esbl_found:
                esbl_pos_patient_data.append(ab_vector)
            else:
                esbl_neg_patient_data.append(ab_vector)

    print("Training on:")
    print("ESBL Positive patients: " + str(len(esbl_pos_patient_data)))
    print("ESBL Negative patients: " + str(len(esbl_neg_patient_data)))
    return esbl_pos_patient_data, esbl_neg_patient_data

def generate_data(patient_data, ab_data, ab_names, AB_CULTURE_COUNT_CUTOFF, CULTURE_SIZE_CUTOFF, ESBL_AB_RESISTENCE_LIST, esbl_result_format=None, numeric=True):
    
    esbl_pos_patient_data, esbl_neg_patient_data = generate_esbl_patient_data(patient_data, ab_data, CULTURE_SIZE_CUTOFF, ESBL_AB_RESISTENCE_LIST)

    RELATIVE_AB_CULTURE_COUNT_CUTOFF = float(len(esbl_pos_patient_data))*(float(AB_CULTURE_COUNT_CUTOFF)/100)

    esbl_pos_ab_count = count_ab(esbl_pos_patient_data, ab_names)
    esbl_neg_ab_count = count_ab(esbl_neg_patient_data, ab_names)
    relevant_ab_list = relevant_ab(esbl_pos_ab_count, ab_data, ab_names, RELATIVE_AB_CULTURE_COUNT_CUTOFF)

    esbl_pos_patient_data = filter_ab(esbl_pos_patient_data, relevant_ab_list, ab_data, esbl_result=esbl_result_format, numeric=numeric)
    esbl_neg_patient_data = filter_ab(esbl_neg_patient_data, relevant_ab_list, ab_data, esbl_result=esbl_result_format, numeric=numeric)

    return esbl_pos_patient_data, esbl_neg_patient_data

def generate_training_data(data_representation, break_amount=1, split_percentage=20):
    while(True):
        train_data = []
        pos_test_data = []
        neg_test_data = []

        if break_amount == 1:
            train_data, pos_test_data = percentage_split([], [],
                                                         data_representation.esbl_pos_patient_data,
                                                         sample_count=int(break_amount*len(data_representation.esbl_pos_patient_data)),
                                                         split_percentage=split_percentage)
        else:
            train_data_temp, pos_test_data = percentage_split([], [],
                                                              data_representation.esbl_pos_patient_data,
                                                              sample_count=int(break_amount*len(data_representation.esbl_pos_patient_data)),
                                                              random_sampling=True, split_percentage=split_percentage)
            train_data = []
            for row in train_data_temp:
                if row not in pos_test_data:
                    train_data.append(row)
            pos_test_data.sort()
            pos_test_data = list(k for k,_ in it.groupby(pos_test_data))

        # testing and training data needs to exist
        if ((not pos_test_data or not train_data) and (split_percentage > 0 and split_percentage < 100)): continue

        # Data exists so append labels for all positive results
        train_labels = ["POS" for i in range(len(train_data))]
        train_data, neg_test_data = percentage_split(train_data, [], data_representation.esbl_neg_patient_data, sample_count=len(train_labels), random_sampling=True, split_percentage=split_percentage)            
 
        # control testing data needs to exist
        if (not neg_test_data and (split_percentage > 0 and split_percentage < 100)): continue

        # Data exists so append leftover labels
        for i in range(len(train_data)-len(train_labels)):
            train_labels.append("NEG")


        # True if data can be used for analysis
        if (train_labels):
            break

    return train_data, pos_test_data, neg_test_data, train_labels

    