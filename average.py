# Libraries
import numpy as np
import pandas as pd
from datetime import datetime
# Modules
import data_reader
import process_data

def init_array(length):
    return 


def run(filename, CULTURE_SIZE_CUTOFF, AB_CULTURE_COUNT_CUTOFF, ESBL_AB_RESISTENCE_LIST):
    csv_data = data_reader.read_csv(filename)
    
    print("Converting data ...")
    id_dict = csv_data[0]
    ab_dict = csv_data[1]
    ab_names = ab_dict.keys()
    print(ab_names)

    esbl_pos_patient_data, esbl_neg_patient_data = process_data.generate_esbl_patient_data(id_dict, ab_dict, CULTURE_SIZE_CUTOFF, ESBL_AB_RESISTENCE_LIST)

    esbl_pos_ab_count = process_data.count_ab(esbl_pos_patient_data, ab_names)
    esbl_neg_ab_count = process_data.count_ab(esbl_neg_patient_data, ab_names)

    print("27", esbl_pos_ab_count)
    print(esbl_pos_ab_count)

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

    relevant_ab_list = process_data.relevant_ab(esbl_pos_ab_count, AB_CULTURE_COUNT_CUTOFF)
    esbl_pos_patient_data = process_data.filter_ab(esbl_pos_patient_data, relevant_ab_list, "POS")
    esbl_neg_patient_data = process_data.filter_ab(esbl_neg_patient_data, relevant_ab_list, "NEG")

    train_data = []
    test_data = []
    train_data, test_data = process_data.percentage_split(train_data, test_data, esbl_pos_patient_data)
    train_data, test_data = process_data.percentage_split(train_data, test_data, esbl_neg_patient_data)

    relevant_ab_names = []
    for ab_index in relevant_ab_list:
        relevant_ab_names.append(ab_names[ab_index])
    relevant_ab_names.append("ESBL")
    print(relevant_ab_names)

    test_dataset = pd.DataFrame(test_data, columns=relevant_ab_names)
    # print(test_dataset.head(200))
    
    # print
    # print("RESULTS OF AVERAGE METHOD:")
    # for index in range(len(esbl_pos_result)):
    #     print("---------------------------")
    #     print("POS", esbl_pos_result[index])
    #     print("NEG", esbl_neg_result[index])