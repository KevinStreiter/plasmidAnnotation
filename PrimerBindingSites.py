import csv
from Bio import SeqIO
from Bio import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqFeature import CompoundLocation
from Bio.SeqFeature import SeqFeature, FeatureLocation
from PrimerRepository import PrimerRepository

class PrimerBindingSites:

    def extractFeatures(self, plasmid_records, primer_binding_sites):

        def writeGeneBankFile(plasmid_records):
            output_file = open('primer_binding_sites.gb', 'w')
            for record in plasmid_records:
                SeqIO.write(record, output_file, 'genbank')
            output_file.close()
            print("---------File:" + " primer_binding_sites.gb " + "created---------")
        
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
        
        def appendFeatures(plasmid_object, start, end, strand, type, qualifier, location_operator):
            if start > end:
                plasmid_object.features.append(SeqFeature(CompoundLocation([FeatureLocation(start, len(plasmid_object.seq)), FeatureLocation(0, end)]),
                                                  type=type, strand=strand,qualifiers=qualifier, location_operator=location_operator))
            else:
                plasmid_object.features.append(SeqFeature(FeatureLocation(start, end), type=type,
                                                  strand=strand,qualifiers=qualifier))

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
                    end_position = evaluateEndPosition(end_position, record, len(record))
                    new_sequence = PrimerRepository(site.id, site.seq, start_position, end_position, forward_strand)
                    list_of_sequence.append(new_sequence)
                    qualifier = {"note": site.id}
                    appendFeatures(record, start_position, end_position, forward_strand, "primer_bind", qualifier, 'join')

                search = plasmid_seq_to_check.reverse_complement().find(site.seq[-15:])
                if search > -1:
                    start_position = search
                    end_position = start_position + 15
                    end_position = evaluateEndPosition(end_position, record, len(record))
                    new_sequence = PrimerRepository(site.id, site.seq, start_position, end_position, backward_strand)
                    list_of_sequence.append(new_sequence)
                    qualifier = {"note": site.id}
                    appendFeatures(record, start_position, end_position, backward_strand, "primer_bind", qualifier, 'join')

        writeCSVFile(list_of_sequence)
        writeGeneBankFile(plasmid_records)
        
        
     







