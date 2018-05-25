import data_reader
import process_data

class Representation():
    """docstring for Representation"""
    def __init__(self):
        self.patient_dict = {}

    def set_culture_parameters(self, CULTURE_SIZE_CUTOFF=0,
                                     AB_CULTURE_COUNT_CUTOFF=0,
                                     ESBL_AB_RESISTENCE_LIST = ["cefotaxim", "ceftazidim", "ceftriaxon"]):
        self.CULTURE_SIZE_CUTOFF = CULTURE_SIZE_CUTOFF
        self.AB_CULTURE_COUNT_CUTOFF = AB_CULTURE_COUNT_CUTOFF
        self.ESBL_AB_RESISTENCE_LIST = ESBL_AB_RESISTENCE_LIST

    def set_RELATIVE_AB_CULTURE_COUNT_CUTOFF(self):
        self.RELATIVE_AB_CULTURE_COUNT_CUTOFF = float(len(self.esbl_pos_patient_data))*(float(self.AB_CULTURE_COUNT_CUTOFF)/100)

    def load_esbl_patient_data(self):
        self.esbl_pos_patient_data, self.esbl_neg_patient_data = process_data.generate_esbl_patient_data(self.id_dict,
                                                                                                         self.ab_dict,
                                                                                                         self.CULTURE_SIZE_CUTOFF,
                                                                                                         self.ESBL_AB_RESISTENCE_LIST)

    def load_filtered_esbl_patient_data(self, esbl_result_format=None, numeric=True):
        self.esbl_pos_patient_data, self.esbl_neg_patient_data = process_data.generate_data(self.id_dict,
                                                                                            self.ab_dict,
                                                                                            self.ab_names,
                                                                                            self.AB_CULTURE_COUNT_CUTOFF,
                                                                                            self.CULTURE_SIZE_CUTOFF,
                                                                                            self.ESBL_AB_RESISTENCE_LIST,
                                                                                            esbl_result_format=esbl_result_format,
                                                                                            numeric=numeric)

    def load_culture_data(self, filename):
        self.csv_data = data_reader.read_csv(filename)
        self.id_dict = self.csv_data[0]
        self.ab_dict = self.csv_data[1]
        self.ab_names = self.ab_dict.keys()