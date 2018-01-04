class SequenceRepository:

    def __init__(self, sequence, feature_type, qualifiers, common_qualifier, qualifier_values):
        self.sequence = sequence
        self.feature_type = feature_type
        self.qualifiers = []
        self.qualifier_values = []
        for qualifier in qualifiers:
            self.qualifiers.append(qualifier)
        self.count = 1
        self.common_qualifier = common_qualifier
        for qualifier_value in qualifier_values:
            self.qualifier_values.append(qualifier_value)

    def getSequence(self):
        return self.sequence

    def getFeature_type(self):
        return self.feature_type

    def getQualifiers(self):
        return self.qualifiers

    def getCount(self):
        return self.count

    def incrementCount(self):
        self.count= self.count + 1

    def appendQualifiers(self,qualifier_list):
        for qualifier in qualifier_list:
            self.qualifiers.append(qualifier)

    def appendQualifierValues(self,qualifier_value_list,key):
        for qualifier_value in qualifier_value_list:
            self.qualifier_values.append(key+":"+qualifier_value)

    def setCommonQualifier(self,qualifier):
        self.common_qualifier = qualifier

    def getQualifierValues(self):
        return self.qualifier_values

    def getCommonQualifier(self):
        return self.common_qualifier




