# -*- coding: utf-8 -*-
"""
@author: Kevin Streiter & Andreas Ott
"""

from Bio import SeqIO
from SequenceRepository import *

vectors = SeqIO.parse("vectors-100.gb", "genbank")



list_of_feature_types = ['promoter', 'oriT', 'rep_origin', 'primer_bind', 'terminator', 'misc_signal',
                         'misc_recomb',
                         'LTR', 'enhancer', '-35_signal', '-10_signal', 'RBS', 'polyA_signal',
                         'sig_peptide', 'CDS',
                         'protein_bind', 'misc_binding', 'mobile_element', 'mRNA', 'tRNA', 'rRNA']

list_of_sequence = []


for record in vectors:
    if len(record.seq) > 1500:
        for feature in record.features:
            for type in list_of_feature_types:
                if feature.type == type:
                    check = False
                    if len(list_of_sequence) > 0:
                        for sequence in list_of_sequence:
                            if sequence.getSequence() == record.seq and sequence.getFeature_type() == feature.type and sequence.getQualifier() == feature.qualifiers.keys():
                                sequence.incrementCount()
                                print "Test"
                                check = True
                    if not check:
                        new_sequence = SequenceRepository(record.seq, feature.type, feature.qualifiers.keys())
                        list_of_sequence.append(new_sequence)








for sequence in list_of_sequence:
    print sequence.getFeature_type()
    print sequence.getSequence()
    print sequence.getCount()
    print sequence.getQualifier()