import data_reader
import process_data
import patient_data
import numpy as np

from copy import deepcopy

class Representation():
    """docstring for Representation"""
    def __init__(self):
        self.patient_dict = {}

    def set_culture_parameters(self, CULTURE_SIZE_CUTOFF=0,
                                     AB_CULTURE_COUNT_CUTOFF=0,
                                     ESBL_AB_RESISTANCE_LIST = ["cefotaxim", "ceftazidim", "ceftriaxon"]):
        self.CULTURE_SIZE_CUTOFF = CULTURE_SIZE_CUTOFF
        self.AB_CULTURE_COUNT_CUTOFF = AB_CULTURE_COUNT_CUTOFF
        self.ESBL_AB_RESISTANCE_LIST = ESBL_AB_RESISTANCE_LIST

    def set_RELATIVE_AB_CULTURE_COUNT_CUTOFF(self):
        self.RELATIVE_AB_CULTURE_COUNT_CUTOFF = float(len(self.esbl_pos_patient_data))*(float(self.AB_CULTURE_COUNT_CUTOFF)/100)

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
                                                                                            numeric=numeric)

    def load_culture_data(self, filename):
        self.id_dict, self.ab_dict = data_reader.read_csv(filename)
        self.ab_names = self.ab_dict.keys()

    def load_patient_data(self, culture_filename, medication_filename):
        # Set parameters
        self.set_culture_parameters()

        # Read culture data
        self.load_culture_data(culture_filename)

        # Process culture data
        self.load_esbl_patient_data()
        self.set_RELATIVE_AB_CULTURE_COUNT_CUTOFF()

        # Extract relevant antibiotics
        relevant_ab_names = process_data.relevant_ab(process_data.count_ab(self.esbl_pos_patient_data, self.ab_names),
                                                                           self.ab_dict,
                                                                           self.ab_names,
                                                                           self.RELATIVE_AB_CULTURE_COUNT_CUTOFF)

        # Read medication data
        medication_dict, temp = data_reader.read_csv(medication_filename, data_type="dot")

        ab_length = len(self.ab_dict.keys())
        for patient_id in self.id_dict.keys():

            if not medication_dict.has_key(patient_id):
                continue

            patient = patient_data.Patient()

            ab_vector, resistance_found, resistance_date = process_data.generate_ab_vector(self.id_dict[patient_id],
                                                                                           ab_length,
                                                                                           self.ESBL_AB_RESISTANCE_LIST,
                                                                                           self.ab_dict,
                                                                                           return_date=True,
                                                                                           numeric=True)

            # Check if data exists
            if np.count_nonzero(ab_vector) < 1:
                continue

            # Set patient data
            patient.patient_id = patient_id
            patient.is_resistant = resistance_found
            patient.resistances = deepcopy(ab_vector)

            medication_vector = np.zeros(ab_length)
            
            # Create DOT vector from ab_vector
            for ab_index in range(len(self.ab_names)):
                ab_name = self.ab_names[ab_index]
                for medication_name in medication_dict[patient_id].keys():
                    if ab_name in medication_name:
                        # DOT is amount of days of medication

                        medication_vector[self.ab_dict[ab_name]] = process_data.count_dot(medication_dict[patient_id][medication_name], limit=resistance_date)

            # Set DOT vector
            patient.medication_use = medication_vector

            self.patient_dict[patient_id] = patient

if __name__ == '__main__':
    x = Representation()
    x.load_patient_data("data/sample2.csv", "data/sample.csv")
    patient_data = x.patient_dict
    random_patient = patient_data.keys()[342]
    print(patient_data[random_patient])