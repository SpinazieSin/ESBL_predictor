import data_reader
import process_data
import patient_data
import numpy as np

from copy import deepcopy

class Representation():
    """docstring for Representation"""
    def __init__(self):
        self.patient_dict = {}
        self.date_range=[2, 3+(6*30)]

    def set_culture_parameters(self, CULTURE_SIZE_CUTOFF = 10,
                                     AB_CULTURE_COUNT_CUTOFF = 50,
                                     ESBL_AB_RESISTANCE_LIST = ["cefotaxim", "ceftazidim", "ceftriaxon", "cefepime"],
                                     RELEVANT_MO_LIST = []):
        self.CULTURE_SIZE_CUTOFF = CULTURE_SIZE_CUTOFF
        self.AB_CULTURE_COUNT_CUTOFF = AB_CULTURE_COUNT_CUTOFF
        self.ESBL_AB_RESISTANCE_LIST = ESBL_AB_RESISTANCE_LIST
        self.RELEVANT_MO_LIST=RELEVANT_MO_LIST

    def set_RELATIVE_AB_CULTURE_COUNT_CUTOFF(self):
        self.RELATIVE_AB_CULTURE_COUNT_CUTOFF = float(len(self.esbl_pos_patient_data))*(float(self.AB_CULTURE_COUNT_CUTOFF)/100)
        # Extract relevant antibiotics
        temp_esbl_pos, temp = process_data.generate_esbl_patient_data(self.id_dict,
                                                                      self.ab_dict,
                                                                      self.CULTURE_SIZE_CUTOFF,
                                                                      self.ESBL_AB_RESISTANCE_LIST,
                                                                      date_range=self.date_range)
        self.relevant_ab_names = process_data.relevant_ab(process_data.count_ab(temp_esbl_pos, self.ab_names),
                                                                                self.ab_dict,
                                                                                self.ab_names,
                                                                                self.RELATIVE_AB_CULTURE_COUNT_CUTOFF)
        # Hard coded relevant ab.
        # self.relevant_ab_names = ["cefotaxim", "ceftazidim", "cefepime"]
        # self.relevant_ab_names = ["tazocin", "ciprofloxacine", "amoxicilline", "cefotaxim", "augmentin", "trimethoprim"]
        self.relevant_ab_names = ["augmentin", "cotrimoxazol", "ciprofloxacine", "trimethoprim", "tazocin", "cefotaxim", "ceftazidim", "cefuroxim", "cefepime", "amoxicilline"]

    def load_esbl_patient_data(self):
        self.esbl_pos_patient_data, self.esbl_neg_patient_data = process_data.generate_esbl_patient_data(self.id_dict,
                                                                                                         self.ab_dict,
                                                                                                         self.CULTURE_SIZE_CUTOFF,
                                                                                                         self.ESBL_AB_RESISTANCE_LIST)

    def load_filtered_esbl_patient_data(self, esbl_result_format=None, numeric=True):
        self.esbl_pos_patient_data, self.esbl_neg_patient_data = process_data.generate_data(self.id_dict,
                                                                                            self.ab_dict,
                                                                                            self.ab_names,
                                                                                            self.AB_CULTURE_COUNT_CUTOFF,
                                                                                            self.CULTURE_SIZE_CUTOFF,
                                                                                            self.ESBL_AB_RESISTANCE_LIST,
                                                                                            esbl_result_format=esbl_result_format,
                                                                                            numeric=numeric,
                                                                                            date_range=self.date_range)
        self.set_RELATIVE_AB_CULTURE_COUNT_CUTOFF()

    def load_culture_data(self, filename):
        self.id_dict, self.ab_dict = data_reader.read_csv(filename, RELEVANT_MO_LIST=self.RELEVANT_MO_LIST)
        self.ab_names = self.ab_dict.keys()

    def load_patient_data(self, culture_filename, medication_filename):
        # Set parameters
        self.set_culture_parameters()

        # Read culture data
        self.load_culture_data(culture_filename)

        # Process culture data
        # self.load_esbl_patient_data()
        # self.set_RELATIVE_AB_CULTURE_COUNT_CUTOFF()
        self.load_filtered_esbl_patient_data()

        # Read medication data
        self.medication_dict, temp = data_reader.read_csv(medication_filename, RELEVANT_MO_LIST=self.RELEVANT_MO_LIST, data_type="dot")

    def generate_patient_data(self):
        ab_length = len(self.ab_dict.keys())
        relevant_ab_length = len(self.relevant_ab_names)

        max_medication_value = np.ones(relevant_ab_length)

        for patient_id in self.id_dict.keys():

            if not self.medication_dict.has_key(patient_id):
                continue

            patient = patient_data.Patient()

            temp_ab_vector, resistance_found, resistance_date = process_data.generate_ab_vector(self.id_dict[patient_id],
                                                                                                   ab_length,
                                                                                                   self.ESBL_AB_RESISTANCE_LIST,
                                                                                                   self.ab_dict,
                                                                                                   return_date=True,
                                                                                                   numeric=True,
                                                                                                   date_range=self.date_range)
            ab_vector = [0 for _ in range(relevant_ab_length)]
            for ab_index in range(relevant_ab_length):
                ab_vector[ab_index] = temp_ab_vector[self.ab_dict[self.relevant_ab_names[ab_index]]]

            # Check if data exists
            if np.count_nonzero(ab_vector) < 1:
                continue

            # Set patient data
            patient.patient_id = patient_id
            patient.is_resistant = resistance_found
            patient.resistances = deepcopy(ab_vector)

            medication_vector = np.zeros(relevant_ab_length)
            
            # Create DOT vector from ab_vector
            for ab_index in range(relevant_ab_length):
                ab_name = self.relevant_ab_names[ab_index]
                for medication_name in self.medication_dict[patient_id].keys():
                    if ab_name in medication_name:
                        # DOT is amount of days of medication
                        dot = process_data.count_dot(self.medication_dict[patient_id][medication_name], limit=resistance_date, date_range=self.date_range)
                        medication_vector[ab_index] = dot
                        # if max_medication_value[ab_index] < dot:
                        #     max_medication_value[ab_index] = dot

            # Set DOT vector
            patient.medication_use = medication_vector

            self.patient_dict[patient_id] = patient

        # SETTING: Normalize dot vectors
        # for patient_id in self.patient_dict.keys():
        #     self.patient_dict[patient_id].medication_use/=max_medication_value

if __name__ == '__main__':
    x = Representation()
    x.set_culture_parameters()
    x.load_culture_data("data/niet_ampc.csv")
    x.load_filtered_esbl_patient_data()    
    print(len(x.id_dict.keys()))
