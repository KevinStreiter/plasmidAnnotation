class SequenceRepository:

    def __init__(self, sequence, feature_type, qualifier):
        self.sequence = sequence
        self.feature_type = feature_type
        self.qualifier = qualifier
        self.count = 1



    def getSequence(self):
        return self.sequence

    def getFeature_type(self):
        return self.feature_type

    def getQualifier(self):
        return self.qualifier

    def getCount(self):
        return self.count

    def incrementCount(self):
        self.count= self.count + 1





