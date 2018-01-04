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
        
        for record in plasmid_records:
            for site in primer_binding_sites:
                search = record.seq.find(site.seq[-15:])
                if search > -1:
                    start_position = search
                    end_position = start_position + 15
                    end_position = evaluateEndPosition(end_position, record, len(record))
                    new_sequence = PrimerRepository(site.id, site.seq, start_position, end_position, forward_strand)
                    list_of_sequence.append(new_sequence)
                    qualifier = {"note": site.id}
                    appendFeatures(record, start_position, end_position, forward_strand, "primer_bind", qualifier, 'join')

                search = record.seq.reverse_complement().find(site.seq[-15:])
                if search > -1:
                    start_position = search
                    end_position = start_position + 15
                    end_position = evaluateEndPosition(end_position, record, len(record))
                    new_sequence = PrimerRepository(site.id, site.seq, start_position, end_position, backward_strand)
                    list_of_sequence.append(new_sequence)
                    qualifier = {"note": site.id}
                    appendFeatures(record, start_position, end_position, backward_strand, "primer_bind", qualifier, 'join')


            firststrand_firstframe = ""
            firststrand_secondframe = ""
            firststrand_thirdframe = ""
            secondstrand_firstframe = ""
            secondstrand_secondframe = ""
            secondstrand_thirdframe = ""
            record.alphabet = IUPAC

            for i in range(0, len(record.seq), 3):
                translated_plasmid = Seq.translate(record.seq[i:i + 3])
                firststrand_firstframe += (str(translated_plasmid))
                translated_reversed_plasmid = Seq.translate(str(record.seq.reverse_complement()[i:i + 3]))
                secondstrand_firstframe += translated_reversed_plasmid
            for i in range(1, len(record.seq), 3):
                translated_plasmid = Seq.translate(record.seq[i:i + 3])
                firststrand_secondframe += (translated_plasmid)
                translated_reversed_plasmid = Seq.translate(str(record.seq.reverse_complement()[i:i + 3]))
                secondstrand_secondframe += translated_reversed_plasmid
            for i in range(2, len(record.seq), 3):
                translated_plasmid = Seq.translate(record.seq[i:i + 3])
                firststrand_thirdframe += (translated_plasmid)
                translated_reversed_plasmid = Seq.translate(str(record.seq.reverse_complement()[i:i + 3]))
                secondstrand_thirdframe += translated_reversed_plasmid

            protein_strands = [Seq.Seq(str(firststrand_firstframe)), Seq.Seq(str(firststrand_secondframe)),
                               Seq.Seq(str(firststrand_thirdframe)), Seq.Seq(str(secondstrand_firstframe)),
                               Seq.Seq(str(secondstrand_secondframe)), Seq.Seq(str(secondstrand_thirdframe))]
            for protein_strand in protein_strands:
                for epitope in special_features:
                    doubled_protein_strand = protein_strand + protein_strand
                    search = doubled_protein_strand.find(epitope.seq)
                    if search > -1:
                        start_position = search*3
                        end_position = start_position + len(epitope.seq) * 3
                        end_position = evaluateEndPosition(end_position, protein_strand, len(protein_strand))
                        print "Found --> ", epitope.id
                        qualifier = {"note": epitope.id}
                        appendFeatures(record, start_position, end_position, 1, "misc_feature", qualifier, 'join')

        writeCSVFile(list_of_sequence)
        writeGeneBankFile(plasmid_records)
        
        
     







