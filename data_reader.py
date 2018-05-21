# Libraries
import csv
import io
from datetime import datetime

def read_csv(filename):
    print("Reading " + filename + " ...")

    with open(filename, 'r') as csvfile:
        
        culture_data = csv.reader(csvfile, delimiter='\t')

        id_dict = {}
        ab_dict = {}

        for culture in culture_data:
            psuedoID = culture[0]
            date = datetime.strptime(culture[1], "%d-%m-%Y")
            ab_name = culture[2]
            if (not ab_dict.has_key(ab_name)):
                ab_dict[ab_name] = len(ab_dict)
            result = culture[5]
            bepaling = culture[8]
            if (not id_dict.has_key(psuedoID)):
                id_dict[psuedoID] = [[date, ab_name, result, bepaling]]
            else:
                id_dict[psuedoID].append([date, ab_name, result, bepaling])
    return id_dict, ab_dict