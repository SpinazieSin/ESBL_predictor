
class Patient():
    def __init__(self):
        self.patient_id = None
        self.is_resistant = None
        self.resistances = None
        self.medication_use = None

    def __str__(self):
        output = ""
        output+="Patient id: {} \n".format(self.patient_id)
        output+="Is resistant: {} \n".format(self.is_resistant)
        output+="Resistance vector: {} \n".format(self.resistances)
        output+="Medication vector: {} \n".format(self.medication_use)
        return output