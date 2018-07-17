# Libraries
import numpy as np
import itertools as it
from random import randint
from math import floor

# Modules
import data
import patient_data

# Important modifications can be found by looking for the SETTING: tag.
# These can be changed by the user to influence the analysis but it is not needed.

def count_ab(esbl_patient_data, ab_names):
    ab_count = np.zeros((len(ab_names), 2))
    for ab_vector in esbl_patient_data:
        if np.count_nonzero(ab_vector) > 0:
            for ab_index in range(len(ab_vector)):
                if ab_vector[ab_index] is 'S':
                    # Suscebtible count
                    ab_count[ab_index][0]+=1
                    # Resistant count
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
                # SETTING: value is the default value of an ab in the ab_vector.
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
            # SETTING: Random probability can be changed to be more/less accurate
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

def find_esbl_pos_day(patient_data, ESBL_AB_RESISTANCE_LIST):
    for culture in patient_data:
        ab = culture[1]
        result = culture[2]
        bepaling = culture[3]
        # SETTING: Culture type can be changed from blood culture (M_banaal_BK) to other types.
        # This only affects the check if a patient is infected with a resistant virus.
        if (ab in ESBL_AB_RESISTANCE_LIST and result is not 'S' and "M_banaal_BK" in bepaling):
            return culture[0]
    return None

def generate_ab_vector(patient_data, ab_length, ESBL_AB_RESISTANCE_LIST, ab_dict, return_date=False, numeric=False, date_range=[5, 90]):
    esbl_found = False
    # SETTING: Non-reversed makes the analyzer take the most recent culture for analysis (which is generally prefered). Reversing
    # it makes it take the oldest culture.
    patient_data.sort(reverse=False)

    if numeric:
        ab_vector = np.zeros(ab_length)
    else:
        ab_vector = [None for x in range(ab_length)]

    # SETTING: Patients with little data can be skipped by increasing the required length of the patient data.
    # This is done later as well with CULTURE_SIZE_CUTOFF, but the minimum length is more basic and saves processing time. 
    if len(patient_data) < 1: return

    esbl_pos_day = find_esbl_pos_day(patient_data, ESBL_AB_RESISTANCE_LIST)
    if esbl_pos_day:
        # The day the esbl is found is the cutoff day
        cutoff_day = esbl_pos_day
    else:
        # Last day day is the cutoff day
        cutoff_day = patient_data[len(patient_data)-1][0] 

    for culture in patient_data:

        culture_day = culture[0]
        ab = culture[1]
        result = culture[2]

        if culture_day == esbl_pos_day:
            esbl_found = True

        total_days = (cutoff_day-culture_day).days
        if total_days > date_range[1]:
            continue

        if esbl_pos_day is not None:
            if total_days < date_range[0]:
                continue

        if result is 'S':
            if numeric:
                ab_vector[ab_dict[ab]] = 1
            else:
                ab_vector[ab_dict[ab]] = 'S'
        else:
            if numeric:
                ab_vector[ab_dict[ab]] = -1
            else:
                ab_vector[ab_dict[ab]] = 'R'
    if return_date:
        return ab_vector, esbl_found, esbl_pos_day
    else:
        return ab_vector, esbl_found

def generate_esbl_patient_data(id_dict, ab_dict, CULTURE_SIZE_CUTOFF, ESBL_AB_RESISTANCE_LIST, date_range=[5, 90]):
    esbl_pos_patient_data = []
    esbl_neg_patient_data = []
    ab_length = len(ab_dict.keys())
    
    for patient in id_dict.keys():

        patient_data = id_dict[patient]

        ab_vector, esbl_found = generate_ab_vector(patient_data, ab_length, ESBL_AB_RESISTANCE_LIST, ab_dict, date_range=date_range)

        # PREVIOUS VERSION OF CUTOFF
        # if len(patient_data) < CULTURE_SIZE_CUTOFF:
        #     continue
        
        # NEW VERSION OF CUTOFF
        if np.count_nonzero(ab_vector) > CULTURE_SIZE_CUTOFF:
            # SETTING: r_length sets the minimum amount of ab resistence needed for a patient to be added.
            r_length = len([1 for x in ab_vector if x == 'R'])
            if esbl_found:
                if r_length < 3:
                    # continue
                    pass
                esbl_pos_patient_data.append(ab_vector)
            else:
                esbl_neg_patient_data.append(ab_vector)

    return esbl_pos_patient_data, esbl_neg_patient_data

def generate_data(patient_data, ab_data, ab_names, AB_CULTURE_COUNT_CUTOFF, CULTURE_SIZE_CUTOFF, ESBL_AB_RESISTANCE_LIST, esbl_result_format=None, numeric=True, date_range=[5, 90]):
    
    esbl_pos_patient_data, esbl_neg_patient_data = generate_esbl_patient_data(patient_data, ab_data, CULTURE_SIZE_CUTOFF, ESBL_AB_RESISTANCE_LIST, date_range=date_range)

    RELATIVE_AB_CULTURE_COUNT_CUTOFF = float(len(esbl_pos_patient_data))*(float(AB_CULTURE_COUNT_CUTOFF)/100)

    esbl_pos_ab_count = count_ab(esbl_pos_patient_data, ab_names)
    esbl_neg_ab_count = count_ab(esbl_neg_patient_data, ab_names)

    # Base relevant ab on RELATIVE_AB_CULTURE_COUNT_CUTOFF.
    relevant_ab_list = relevant_ab(esbl_pos_ab_count, ab_data, ab_names, RELATIVE_AB_CULTURE_COUNT_CUTOFF)
    # Manually add relevant ab.
    # relevant_ab_list = ["cefotaxim", "ceftazidim", "cefepime"]
    # relevant_ab_list = ["tazocin", "ciprofloxacine", "amoxicilline", "fosfomycine", "augmentin", "trimethoprim"]
    relevant_ab_list = ["augmentin", "cotrimoxazol", "ciprofloxacine", "trimethoprim", "tazocin", "cefotaxim", "ceftazidim", "cefuroxim", "cefepime", "amoxicilline"]

    esbl_pos_patient_data = filter_ab(esbl_pos_patient_data, relevant_ab_list, ab_data, esbl_result=esbl_result_format, numeric=numeric)
    esbl_neg_patient_data = filter_ab(esbl_neg_patient_data, relevant_ab_list, ab_data, esbl_result=esbl_result_format, numeric=numeric)


    print("Found:")
    print("Positive resistance patients: " + str(len(esbl_pos_patient_data)))
    print("Negative resistance patients: " + str(len(esbl_neg_patient_data)))
    return esbl_pos_patient_data, esbl_neg_patient_data

def generate_training_data(pos_patient_data, neg_patient_data, break_amount=1, split_percentage=20):
    while(True):
        train_data = []
        pos_test_data = []
        neg_test_data = []

        # HACK: Do not mess with this unless it is carefully considered.
        if break_amount == 1:
            train_data, pos_test_data = percentage_split([], [],
                                                         pos_patient_data,
                                                         sample_count=int(break_amount*len(pos_patient_data)),
                                                         split_percentage=split_percentage)
        else:
            train_data_temp, pos_test_data_temp = percentage_split([], [],
                                                              pos_patient_data,
                                                              sample_count=int(break_amount*len(pos_patient_data)),
                                                              random_sampling=True, split_percentage=split_percentage)
            train_data = []
            row_removed = [] # HACK: Added so no extra data is removed.
            pos_test_data = [ list(item) for item in pos_test_data_temp ]
            # TODO: Make rows indexed so we don't need to find them again and use the row_removed hack.
            for temp_row in train_data_temp:
                row = list(temp_row)
                # If train data is in test data, don't use it.
                if row not in pos_test_data or row in row_removed:
                    train_data.append(row)
                # In case of duplicate data.
                else:
                    row_removed.append(row)
            pos_test_data.sort()
            pos_test_data = list(k for k,_ in it.groupby(pos_test_data))

        # testing and training data needs to exist
        if ((not pos_test_data or not train_data) and (split_percentage > 0 and split_percentage < 100)): continue

        # Data exists so append labels for all positive results
        train_labels = ["POS" for i in range(len(train_data))]
        train_data, neg_test_data = percentage_split(train_data, [], neg_patient_data, sample_count=len(train_labels), random_sampling=True, split_percentage=split_percentage)

        # control testing data needs to exist
        if (not neg_test_data and (split_percentage > 0 and split_percentage < 100)): continue

        # Data exists so append leftover labels
        for i in range(len(train_data)-len(train_labels)):
            train_labels.append("NEG")


        # True if data can be used for analysis
        if (train_labels):
            break

    return train_data, pos_test_data, neg_test_data, train_labels

def count_dot(dates, limit=False, date_range=[2, 3+(6*30)]):

    dates.sort(reverse=True)

    if limit:
        init_date = limit
    else:
        init_date = dates[0]

    dot = 0
    for date in dates:
        if date > init_date: continue
        day_diff = (init_date-date).days
        if day_diff > date_range[0] and day_diff < date_range[1]:
            dot+=1

    return dot

def vectorize_medication_data(patient_dict):
    pos_data = []
    neg_data = []

    for patient_id in patient_dict.keys():
        current_patient = patient_dict[patient_id]
        if current_patient.is_resistant:
            pos_data.append(current_patient.medication_use)
        else:
            neg_data.append(current_patient.medication_use)
    return pos_data, neg_data

def vectorize_culture_data(patient_dict):
    pos_data = []
    neg_data = []

    for patient_id in patient_dict.keys():
        current_patient = patient_dict[patient_id]
        if current_patient.is_resistant:
            pos_data.append(current_patient.resistances)
        else:
            neg_data.append(current_patient.resistances)
    return pos_data, neg_data