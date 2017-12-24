"""
@author: Kevin Streiter & Andreas Ott
"""

from SequenceRepository import *
import csv
from collections import Counter


class CommonFeatures:

    def extractFeatures(self, plasmid_records):
        
        def writeCSVFile(list_of_sequence):
            file = csv.writer(open("common_features.csv", "w+"))
            header_keys = ['seq', 'feature_type', 'common_qualifier', 'count']
            file.writerow(header_keys)
            for sequence in list_of_sequence:
                if sequence.getCount() > 2 and len(sequence.getSequence()) >= 4:
                    write_row = list()
                    write_row.append(sequence.getSequence())
                    write_row.append(sequence.getFeature_type())
                    write_row.append(sequence.getCommonQualifier())
                    write_row.append(sequence.getCount())
                    file.writerow(write_row)
                    #file.close()
            print("---------File:" + " common_features.csv " + "created---------")

        list_of_sequence = []

        list_of_feature_types = ['promoter', 'oriT', 'rep_origin', 'primer_bind', 'terminator', 'misc_signal','misc_recomb','LTR', 'enhancer', '-35_signal', '-10_signal', 'RBS', 'polyA_signal', 'sig_peptide', 'CDS','protein_bind', 'misc_binding', 'mobile_element', 'mRNA', 'tRNA', 'rRNA']

        for record in plasmid_records:
            if len(record.seq) > 1500:
                for feature in record.features:
                    for type in list_of_feature_types:
                        if feature.type == type:
                            check = False
                            if len(list_of_sequence) > 0:
                                for sequence in list_of_sequence:
                                    if sequence.getSequence() == feature.extract(record.seq) and sequence.getFeature_type() == feature.type:
                                        sequence.incrementCount()
                                        for qualifier in feature.qualifiers.keys():
                                                sequence.appendQualifierValues(feature.qualifiers[qualifier])
                                        sequence.appendQualifiers(feature.qualifiers.keys())
                                        check = True
                            if not check:
                                qualifier_values = []
                                for qualifier in feature.qualifiers.keys():
                                    qualifier_values.append(qualifier)
                                new_sequence = SequenceRepository(feature.extract(record.seq), feature.type, feature.qualifiers.keys(), '',qualifier_values)
                                #print(new_sequence.getSequence(), new_sequence.getQualifiers(), new_sequence.getFeature_type())
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

        for sequence in list_of_sequence:
            if sequence.getCount() > 2:
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
               
                common_qualifier,num_most_qualifier = Counter(sequence.getQualifierValues()).most_common(1)[0]
                sequence.setCommonQualifier(common_qualifier)
            
        writeCSVFile(list_of_sequence)
