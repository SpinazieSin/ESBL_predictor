# Libraries
import numpy as np
from scipy import stats
from datetime import datetime

# Modules
import data_reader


def find_esbl_pos_day(patient_data, ESBL_AB_RESISTENCE_LIST):
    for culture in patient_data:
        ab = culture[1]
        result = culture[2]
        if (ab in ESBL_AB_RESISTENCE_LIST and result is not 'S'):
            return culture[0]
    return None

def run(filename, CULTURE_SIZE_CUTOFF, AB_CULTURE_COUNT_CUTOFF, ESBL_AB_RESISTENCE_LIST):
    csv_data = data_reader.read_csv(filename)
    
    print("Converting data ...")
    id_dict = csv_data[0]
    ab_dict = csv_data[1]

    ab_names = ab_dict.keys()
    esbl_pos_ab_count = np.zeros((len(ab_names), 2))
    esbl_neg_ab_count = np.zeros((len(ab_names), 2))
    esbl_pos_patient_data = []
    esbl_neg_patient_data = []

    for patient in id_dict.keys():

        esbl_found = False
        patient_data = id_dict[patient]
        if len(patient_data) < CULTURE_SIZE_CUTOFF:
            continue
        
        patient_data.sort()
        ab_vector = np.zeros(len(ab_names))

        esbl_pos_day = find_esbl_pos_day(patient_data, ESBL_AB_RESISTENCE_LIST)

        for culture in patient_data:
            culture_day = culture[0]
            ab = culture[1]
            result = culture[2]
            if culture_day == esbl_pos_day:
                esbl_found = True
                break
            else:
                if result is 'S':
                    ab_vector[ab_dict[ab]] = 1
                else:
                    ab_vector[ab_dict[ab]] = -1


        if esbl_found:
            if np.count_nonzero(ab_vector) > 0:
                esbl_pos_patient_data.append(ab_vector)
                for ab_index in range(len(ab_vector)):
                    ab_result = ab_vector[ab_index]
                    if ab_result > 0:
                        # Suscebtible count
                        esbl_pos_ab_count[ab_index][0]+=1
                        # Resistent count
                        esbl_pos_ab_count[ab_index][1]+=1
                    elif ab_result < 0:
                        # Reistent count
                        esbl_pos_ab_count[ab_index][1]+=1
        else:
            esbl_neg_patient_data.append(ab_vector)
            for ab_index in range(len(ab_vector)):
                ab_result = ab_vector[ab_index]
                if ab_result > 0:
                    # Suscebtible count
                    esbl_neg_ab_count[ab_index][0]+=1
                    # Resistent count
                    esbl_neg_ab_count[ab_index][1]+=1
                elif ab_result < 0:
                    # Reistent count
                    esbl_neg_ab_count[ab_index][1]+=1
        print("----------------------------------------------------")
        print(esbl_pos_day)
        print(ab_vector)
        print(esbl_found)
        print(len(patient_data))


    print("Analysing ...")
    esbl_pos_result = []
    esbl_neg_result = []
    for ab in ab_names:
        ab_index = ab_dict[ab]
        S = esbl_pos_ab_count[ab_index][0]
        SIR = esbl_pos_ab_count[ab_index][1]

        if (SIR > AB_CULTURE_COUNT_CUTOFF):
            esbl_pos_result.append([ab, S/SIR, S, SIR])
            S = esbl_neg_ab_count[ab_index][0]
            SIR = esbl_neg_ab_count[ab_index][1]

            if (SIR > AB_CULTURE_COUNT_CUTOFF):
                esbl_neg_result.append([ab, S/SIR, S, SIR])

    print
    print("RESULTS OF AVERAGE METHOD:")
    for index in range(len(esbl_pos_result)):
        print("---------------------------")
        print("POS", esbl_pos_result[index])
        print("NEG", esbl_neg_result[index])