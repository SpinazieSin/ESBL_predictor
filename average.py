# Libraries
import numpy as np
from datetime import datetime

# Modules
import data
import data_reader
import process_data

class Average():
    def __init__(self, filename, CULTURE_SIZE_CUTOFF=0, AB_CULTURE_COUNT_CUTOFF=0, ESBL_AB_RESISTENCE_LIST=[]):
        self.data = data.Representation(filename)
        self.data.set_culture_parameters(CULTURE_SIZE_CUTOFF=CULTURE_SIZE_CUTOFF, AB_CULTURE_COUNT_CUTOFF=AB_CULTURE_COUNT_CUTOFF)
        self.esbl_pos_ab_count = []
        self.esbl_neg_ab_count = []


    def run(self):
        print("Converting data ...")
        self.data.load_esbl_patient_data()

        self.esbl_pos_ab_count = process_data.count_ab(self.data.esbl_pos_patient_data, self.data.ab_names)
        self.esbl_neg_ab_count = process_data.count_ab(self.data.esbl_neg_patient_data, self.data.ab_names)

        print("Analysing")
        self.analyse()

    def analyse(self):
        print("Analysing ...")
        esbl_pos_result = []
        esbl_neg_result = []
        for ab in self.data.ab_names:
            ab_index = self.data.ab_dict[ab]
            S = self.esbl_pos_ab_count[ab_index][0]
            SIR = self.esbl_pos_ab_count[ab_index][1]

            if (SIR > self.data.RELATIVE_AB_CULTURE_COUNT_CUTOFF):
                esbl_pos_result.append([ab, S/SIR, S, SIR])
                S = self.esbl_neg_ab_count[ab_index][0]
                SIR = self.esbl_neg_ab_count[ab_index][1]

                if (SIR > self.data.RELATIVE_AB_CULTURE_COUNT_CUTOFF):
                    esbl_neg_result.append([ab, S/SIR, S, SIR])
        
        print("RESULTS OF AVERAGE METHOD:")
        for index in range(len(esbl_pos_result)):
            print("________")
            print("Positive resistence result: Name - {}, S/SIR - {}, Susceptible Count - {}, Resistent Count - {}".format(esbl_pos_result[index][0], esbl_pos_result[index][1], esbl_pos_result[index][2], esbl_pos_result[index][3]))
            print("Negative resistence result: Name - {}, S/SIR - {}, Susceptible Count - {}, Resistent Count - {}".format(esbl_neg_result[index][0], esbl_neg_result[index][1], esbl_neg_result[index][2], esbl_neg_result[index][3]))
