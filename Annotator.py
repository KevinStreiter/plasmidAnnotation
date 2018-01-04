"""
@author: Kevin Streiter & Andreas Ott
"""

from Bio.SeqFeature import CompoundLocation
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO

class Annotator:
        @staticmethod
        def writeGeneBankFile(plasmid_records,filename):
            output_file = open("OutputFiles/"+filename, 'w')
            for record in plasmid_records:
                SeqIO.write(record, output_file, 'genbank')
            output_file.close()
            print("---------File:" + filename + "created---------")
        @staticmethod
        def appendFeatures(plasmid_object, start, end, strand, type, qualifier, location_operator):
            if start > end:
                plasmid_object.features.append(SeqFeature(CompoundLocation([FeatureLocation(start, len(plasmid_object.seq)), FeatureLocation(0, end)]),
                                                          type=type, strand=strand,qualifiers=qualifier, location_operator=location_operator))
            else:
                plasmid_object.features.append(SeqFeature(FeatureLocation(start, end), type=type,
                                                          strand=strand,qualifiers=qualifier))
        @staticmethod
        def evaluateEndPosition(end_position, record, record_length):
            if end_position > record_length - 1:
                end_position = end_position - record_length
                return end_position
            return end_position