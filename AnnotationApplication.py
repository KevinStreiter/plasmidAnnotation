"""
@author: Kevin Streiter & Andreas Ott
"""
import sys
from CommonFeatures import CommonFeatures
from PrimerBindingSites import PrimerBindingSites
from SpecialFeatures import SpecialFeatures
from Bio.Alphabet import IUPAC
from Bio import SeqRecord
from Blaster import *


def main(argv):

    checkSuffixFasta = ".fasta"
    checkSuffixGb = ".gb"
    plasmid_records = []

    if len(argv) == 0:
        message = "Please provide a file"
        sys.exit(message)

    filename = str(argv[0])

    if filename.endswith(checkSuffixFasta):
        input_fasta = SeqIO.parse(filename, "fasta")
        for fastRec in input_fasta:
            fastRec = Seq.Seq(str(fastRec.seq), IUPAC.ambiguous_dna)
            plasmid_records.append(SeqRecord.SeqRecord(fastRec))

    elif filename.endswith(checkSuffixGb):
        plasmid_records = list(SeqIO.parse(filename, "genbank"))

    LEARNING_FILE = "InputFiles/vectors.gb"
    vectors = SeqIO.parse(LEARNING_FILE,"genbank")
    common_primer_file = "common_primer.mfasta"
    special_features = "tags_epitopes.mfasta"

    common_primer_records = list(SeqIO.parse(common_primer_file, "fasta"))
    special_features_records = list(SeqIO.parse(special_features, "fasta"))
    
    sequenceRepository = CommonFeatures().extractFeatures(vectors)
    CommonFeatures().annotateSequences(plasmid_records,sequenceRepository)
    PrimerBindingSites().extractFeatures(plasmid_records, common_primer_records)
    SpecialFeatures().extractFeatures(plasmid_records, special_features_records)
    Blaster().blastSearch(plasmid_records)
    
if __name__ == "__main__":
    main(sys.argv[1:])












