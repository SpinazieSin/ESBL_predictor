# Libraries
import numpy as np
from random import randint

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

def relevant_ab(esbl_pos_ab_count, ab_dict, ab_names, AB_CULTURE_COUNT_CUTOFF):
    relevant_ab_list = []
    for ab in ab_names:
        ab_index = ab_dict[ab]
        if esbl_pos_ab_count[ab_index][1] > AB_CULTURE_COUNT_CUTOFF:
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

def percentage_split(train_data, test_data, esbl_patient_data, split_percentage=10):
    for data_row in esbl_patient_data:
        if randint(0, 100) <= split_percentage:
            test_data.append(data_row)
        else:
            train_data.append(data_row)
    return train_data, test_data 

def find_esbl_pos_day(patient_data, ESBL_AB_RESISTENCE_LIST):
    for culture in patient_data:
        ab = culture[1]
        result = culture[2]
        if (ab in ESBL_AB_RESISTENCE_LIST and result is not 'S'):
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

        for culture in patient_data:
            culture_day = culture[0]
            ab = culture[1]
            result = culture[2]
            if culture_day == esbl_pos_day:
                esbl_found = True
                # Remove this break if you want to include the cultures of the day
                # where the patient had the ESBL producing microorganism.
                break
            else:
                if result is 'S':
                    ab_vector[ab_dict[ab]] = 'S'
                else:
                    ab_vector[ab_dict[ab]] = 'R'

        if np.count_nonzero(ab_vector) > 0:
            if esbl_found:
                esbl_pos_patient_data.append(ab_vector)
            else:
                esbl_neg_patient_data.append(ab_vector)
    return esbl_pos_patient_data, esbl_neg_patient_data