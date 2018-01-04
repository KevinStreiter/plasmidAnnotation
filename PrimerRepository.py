class PrimerRepository:
    
    def __init__(self, primer_name, sequence, start_position, end_position, strand):
        self.primer_name = primer_name
        self.sequence = sequence
        self.start_position = start_position
        self.end_position = end_position
        self.strand = strand

    def getPrimerName(self):
        return self.primer_name
    
    def getSequence(self):
        return self.sequence

    def getStartPosition(self):
        return self.start_position

    def getEndPosition(self):
        return self.end_position

    def getStrand(self):
        return self.strand


    

    
    
    