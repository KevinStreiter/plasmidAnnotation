import csv
from Bio import SeqIO
from Bio import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqFeature import CompoundLocation
from Bio.SeqFeature import SeqFeature, FeatureLocation
from PrimerRepository import PrimerRepository
from Annotator import *

class PrimerBindingSites:

    def extractFeatures(self, plasmid_records, primer_binding_sites):

        def writeCSVFile(list_of_sequence):
            file = csv.writer(open("primer_binding_sites.csv", "w+"))
            header_keys = ['primer_id', 'primer_seq', 'start_position', 'end_position','strand']
            file.writerow(header_keys)
            for sequence in list_of_sequence:
                write_row = list()
                write_row.append(sequence.getPrimerName())
                write_row.append(sequence.getSequence())
                write_row.append(sequence.getStartPosition())
                write_row.append(sequence.getEndPosition())
                write_row.append(sequence.getStrand())
                file.writerow(write_row)
                #file.close()
            print("---------File:" + " primer_binding_sites.csv " + "created---------")


        forward_strand = 1
        backward_strand = -1
        
        list_of_sequence = []
        longest_primer_site = 0
        for site in primer_binding_sites:
            if len(site.seq)>longest_primer_site:
                longest_primer_site = len(site.seq)

        for record in plasmid_records:
            #Seq appended to check for Sites over "0"
            plasmid_seq_to_check = record.seq + record.seq[:longest_primer_site]
            for site in primer_binding_sites:
                search = plasmid_seq_to_check.find(site.seq[-15:])
                if search > -1:
                    start_position = search
                    end_position = start_position + 15
                    end_position = Annotator().evaluateEndPosition(end_position, record, len(record))
                    new_sequence = PrimerRepository(site.id, site.seq, start_position, end_position, forward_strand)
                    list_of_sequence.append(new_sequence)
                    qualifier = {"note": site.id}
                    Annotator().appendFeatures(record, start_position, end_position, forward_strand, "primer_bind", qualifier, 'join')

                search = plasmid_seq_to_check.reverse_complement().find(site.seq[-15:])
                if search > -1:
                    start_position = search
                    end_position = start_position + 15
                    end_position = Annotator().evaluateEndPosition(end_position, record, len(record))
                    new_sequence = PrimerRepository(site.id, site.seq, start_position, end_position, backward_strand)
                    list_of_sequence.append(new_sequence)
                    qualifier = {"note": site.id}
                    Annotator().appendFeatures(record, start_position, end_position, backward_strand, "primer_bind", qualifier, 'join')

        writeCSVFile(list_of_sequence)
        Annotator().writeGeneBankFile(plasmid_records,"primer_binding_sites.gb")
        
        
     







