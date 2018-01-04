
from Annotator import *

class PrimerBindingSites:

    def extractFeatures(self, plasmid_records, primer_binding_sites):


        forward_strand = 1
        backward_strand = -1

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
                    qualifier = {"note": site.id}
                    Annotator().appendFeatures(record, start_position, end_position, forward_strand, "primer_bind", qualifier, 'join')

                search = plasmid_seq_to_check.reverse_complement().find(site.seq[-15:])
                if search > -1:
                    start_position = search
                    end_position = start_position + 15
                    end_position = Annotator().evaluateEndPosition(end_position, record, len(record))
                    qualifier = {"note": site.id}
                    Annotator().appendFeatures(record, start_position, end_position, backward_strand, "primer_bind", qualifier, 'join')

        Annotator().writeGeneBankFile(plasmid_records,"primer_binding_sites.gb")
        
        
     







