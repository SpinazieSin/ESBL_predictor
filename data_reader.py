# Libraries
import csv
import io
from datetime import datetime

def read_csv(filename, data_type="culture"):
    print("Reading " + filename + " ...")

    with open(filename, 'r') as csvfile:
        
        csv_data = csv.reader(csvfile, delimiter='\t')

        id_dict = {}
        ab_dict = {}

        if data_type == "culture":

            for culture in csv_data:
                pseudoID = culture[0]
                date = datetime.strptime(culture[1], "%d-%m-%Y").date()
                ab_name = culture[2].lower()
                if (not ab_dict.has_key(ab_name)):
                    ab_dict[ab_name] = len(ab_dict)
                result = culture[5]
                bepaling = culture[8]
                if (not id_dict.has_key(pseudoID)):
                    id_dict[pseudoID] = [[date, ab_name, result, bepaling]]
                else:
                    id_dict[pseudoID].append([date, ab_name, result, bepaling])
        
        elif data_type == "dot":

            for medication in csv_data:
                pseudoID = medication[0]
                date = datetime.strptime(medication[1], "%d-%m-%Y %H:%M:%S").date()
                ab_name = medication[3].lower()
                if (not id_dict.has_key(pseudoID)):
                    id_dict[pseudoID] = {}
                    id_dict[pseudoID][ab_name] = [date]
                else:
                    if (not id_dict[pseudoID].has_key(ab_name)):
                        id_dict[pseudoID][ab_name] = [date]
                    elif (date not in id_dict[pseudoID][ab_name]):
                        id_dict[pseudoID][ab_name].append(date)

    return id_dict, ab_dict

if __name__ == '__main__':
    filename = "data/sample.csv"
    x, y = read_csv(filename, data_type="dot")
    print("Patient")
    print(x.keys()[0])
    print("Medications names")
    print(x["212098"].keys())
    print("Days of the ")
    print(x[x.keys()[0]][x[x.keys()[0]].keys()[1]])