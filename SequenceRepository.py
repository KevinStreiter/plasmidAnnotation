class SequenceRepository:

    def __init__(self, sequence, feature_type, qualifiers):
        self.sequence = sequence
        self.feature_type = feature_type
        self.qualifiers = []
        for qualifier in qualifiers:
            self.qualifiers.append(qualifier)
        self.count = 1



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





