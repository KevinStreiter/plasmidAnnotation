# -*- coding: utf-8 -*-
"""
@author: Kevin Streiter & Andreas Ott
"""

from Bio import SeqIO
from SequenceRepository import *

vectors = SeqIO.parse("vectors-100.gb", "genbank")
list_of_sequence = []

list_of_feature_types = ['promoter', 'oriT', 'rep_origin', 'primer_bind', 'terminator', 'misc_signal','misc_recomb','LTR', 'enhancer', '-35_signal', '-10_signal', 'RBS', 'polyA_signal', 'sig_peptide', 'CDS','protein_bind', 'misc_binding', 'mobile_element', 'mRNA', 'tRNA', 'rRNA']


for record in vectors:
    if len(record.seq) > 1500:
        for feature in record.features:
            for type in list_of_feature_types:
                if feature.type == type:
                    check = False
                    if len(list_of_sequence) > 0:
                        for sequence in list_of_sequence:
                            if sequence.getSequence() == record.seq and sequence.getFeature_type() == feature.type:
                                sequence.incrementCount()
                                sequence.appendQualifiers(feature.qualifiers.keys())
                                check = True
                    if not check:
                        new_sequence = SequenceRepository(record.seq, feature.type, feature.qualifiers.keys())
                        list_of_sequence.append(new_sequence)

note_type = ['promoter', 'oriT', 'rep_origin', 'primer_bind', 'terminator', 'misc_signal','misc_recomb','LTR', 'enhancer', '-35_signal', '-10_signal', 'RBS', 'polyA_signal','sig_peptide','protein_bind', 'misc_binding','mobile_element']
note_list = []
product_type = ['CDS','tRNA', 'rRNA']
product_list  = []
bound_moiety_type = ['protein_bind', 'misc_binding']
bound_moiety_list = []
mobile_element_type = ['mobile_element']
mobile_element_list = []
gene_type = ['mRNA','-35_signal', '-10_signal', 'RBS', 'polyA_signal','sig_peptide', 'CDS']
gene_list = []
undefined_qualifier_list = []
count = 0
for sequence in list_of_sequence:
    #if sequence.getCount() > 2:
    if sequence.getFeature_type() in note_type:
        note_type_bool = False
        for qualifier in sequence.getQualifiers():
            if qualifier == 'note':
                note_type_bool = True
        if note_type_bool:
            note_list.append(sequence)
        else:
            undefined_qualifier_list.append(sequence)
    elif sequence.getFeature_type() in product_type:
        product_type_bool = False
        for qualifier in sequence.getQualifiers():
            if qualifier == 'product':
                product_type_bool = True
        if product_type_bool:
            product_list.append(sequence)
        else:
            undefined_qualifier_list.append(sequence)
    elif sequence.getFeature_type() in bound_moiety_type:
        bound_moiety_bool = False
        for qualifier in sequence.getQualifiers():
            if qualifier == 'bound_moiety':
                bound_moiety_bool = True
        if bound_moiety_bool:
            bound_moiety_list.append(sequence)
        else:
            undefined_qualifier_list.append(sequence)
    elif sequence.getFeature_type() in mobile_element_type:
        mobile_bool = False
        for qualifier in sequence.getQualifiers():
            if qualifier == 'mobile':
                mobile_bool = True
        if mobile_bool:
            mobile_element_list.append(sequence)
        else:
            undefined_qualifier_list.append(sequence)
    elif sequence.getFeature_type() in gene_type:
        gene_bool = False
        for qualifier in sequence.getQualifiers():
            if qualifier == 'gene':
                gene_bool = True
        if gene_bool:
            gene_list.append(sequence)
        else:
            undefined_qualifier_list.append(sequence)

print("Note-List")
for sequence in mobile_element_list:
    print(sequence.getFeature_type())
print (count)
