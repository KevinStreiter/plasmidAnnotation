import csv
from Bio.SeqFeature import CompoundLocation
from Bio.SeqFeature import SeqFeature, FeatureLocation
from PrimerRepository import PrimerRepository

class PrimerBindingSites:

    def extractFeatures(self, plasmid_records, primer_binding_sites):

#        def writeGeneBankFile(plasmid_records):
#            output_file = open('primer_binding_sites.gb', 'w')
#            for record in plasmid_records:
#                SeqIO.write(record, output_file, 'genbank')
#            output_file.close()
#            print("---------File:" + " primer_binding_sites.gb " + "created---------")
        
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
            
        def evaluateEndPosition(end_position, record, record_length):
            if end_position > record_length - 1:
                end_position = end_position - record_length
                return end_position
            return end_position

        forward_strand = 1
        backward_strand = -1
        
        list_of_sequence = []
        
        for record in plasmid_records:
            for site in primer_binding_sites:
                search = record.seq.find(site.seq[-15:])
                if search > -1:
                    start_position = search
                    end_position = start_position + 15
                    end_position = evaluateEndPosition(end_position, record, len(record))
                    new_sequence = PrimerRepository(site.id, site.seq, start_position, end_position, forward_strand)
                    list_of_sequence.append(new_sequence)

                search = record.seq.reverse_complement().find(site.seq[-15:])
                if search > -1:
                    start_position = search
                    end_position = start_position + 15
                    end_position = evaluateEndPosition(end_position, record, len(record))
                    new_sequence = PrimerRepository(site.id, site.seq, start_position, end_position, backward_strand)
                    list_of_sequence.append(new_sequence)

        writeCSVFile(list_of_sequence)
     







