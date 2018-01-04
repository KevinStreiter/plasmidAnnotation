from Bio import Seq
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Annotator import *
from Bio.SeqFeature import CompoundLocation
from Bio.SeqFeature import SeqFeature, FeatureLocation

class SpecialFeatures:
    
    def extractFeatures(self, records, special_features):

        for record in records:

            strand_first = ""
            reverse_strand_first = ""
            strand_second = ""
            reverse_strand_second = ""
            strand_third = ""
            reverse_strand_third = ""
            record.alphabet = IUPAC
            seq_len = len(record.seq)
            for i in range(0, seq_len, 3):
                if i+3<=seq_len:
                    translated_plasmid = Seq.translate(record.seq[i:i + 3])
                    strand_first += (str(translated_plasmid))
                    translated_reversed_plasmid = Seq.translate(str(record.seq.reverse_complement()[i:i + 3]))
                    reverse_strand_first += translated_reversed_plasmid
                if i+4<=seq_len:
                    translated_plasmid = Seq.translate(record.seq[i+1:i + 4])
                    strand_second += (translated_plasmid)
                    translated_reversed_plasmid = Seq.translate(str(record.seq.reverse_complement()[i+1:i + 4]))
                    reverse_strand_second += translated_reversed_plasmid
                if i+5<=seq_len:
                    translated_plasmid = Seq.translate(record.seq[i+2:i + 5])
                    strand_third += (translated_plasmid)
                    translated_reversed_plasmid = Seq.translate(str(record.seq.reverse_complement()[i+2:i + 5]))
                    reverse_strand_third += translated_reversed_plasmid

            protein_strands = [Seq.Seq(str(strand_first)), Seq.Seq(str(reverse_strand_first)),
                               Seq.Seq(str(strand_second)), Seq.Seq(str(reverse_strand_second)),
                               Seq.Seq(str(strand_third)), Seq.Seq(str(reverse_strand_third))]

            for protein_strand in protein_strands:
                for epitope in special_features:
                    protein_strand_to_check = protein_strand + protein_strand[:3]
                    search = protein_strand_to_check.find(epitope.seq)
                    if search > -1:
                        start_position = search*3
                        end_position = start_position + len(epitope.seq) * 3
                        end_position = Annotator().evaluateEndPosition(end_position, protein_strand, len(protein_strand))
                        print ("Found --> ", epitope.id)
                        qualifier = {"note": epitope.id}
                        Annotator().appendFeatures(record, start_position, end_position, 1, "misc_feature", qualifier, 'join')
                    

        Annotator().writeGeneBankFile(records,"special_translated_features.gb")
