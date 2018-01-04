"""
@author: Kevin Streiter & Andreas Ott
"""

from SequenceRepository import SequenceRepository
import csv
from collections import Counter
from Annotator import *

class CommonFeatures:

    def extractFeatures(self, plasmid_records):
        
        def annotateSequences(list_of_sequence):
            for plasmid in plasmid_records:
                for sequenceObject in list_of_sequence:
                    sequence = sequenceObject.getSequence()
                    doubled_plasmid = plasmid.seq + plasmid.seq
                    search = doubled_plasmid.find(sequence)

                    if search > -1:
                        start_pos = search
                        end_pos = start_pos + len(sequence)
                        end_pos = Annotator().evaluateEndPosition(end_pos, plasmid,len(record))
                        split_qualifier = sequenceObject.getCommonQualifier().split(":")
                        if len(split_qualifier)>1:
                            key = split_qualifier[0]
                            value = split_qualifier[1]
                            qualifier = {key: value}
                            Annotator().appendFeatures(plasmid, start_pos, end_pos, 1, sequenceObject.getFeature_type(), qualifier, 'join')

                    search = doubled_plasmid.reverse_complement().find(sequence)
                    if search > -1:
                        start_pos = search
                        end_pos = start_pos + len(sequence)
                        end_pos = Annotator().evaluateEndPosition(end_pos, plasmid,len(record))
                        split_qualifier = sequenceObject.getCommonQualifier().split(":")
                        if len(split_qualifier)>1:
                            key = split_qualifier[0]
                            value = split_qualifier[1]
                            qualifier = {key: value}
                            Annotator().appendFeatures(plasmid, start_pos, end_pos, -1, sequenceObject.getFeature_type(),qualifier, 'join')
            Annotator().writeGeneBankFile(plasmid_records,"common_features.gb")
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
                                                sequence.appendQualifierValues(feature.qualifiers[qualifier],qualifier)
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
                if(len(sequence.getQualifierValues()) > 0):
                    #for values in sequence.getQualifierValues():
                    #    print(values)
                    common_qualifier,num_most_qualifier = Counter(sequence.getQualifierValues()).most_common(1)[0]
                    sequence.setCommonQualifier(common_qualifier)
            
        annotateSequences(list_of_sequence)
