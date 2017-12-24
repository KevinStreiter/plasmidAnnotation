"""
@author: Kevin Streiter & Andreas Ott
"""
import sys
from CommonFeatures import CommonFeatures
from PrimerBindingSites import PrimerBindingSites
from Bio import SeqIO


def main(argv):
    
    if len(argv) == 0:
        message = "Please provide a file"
        sys.exit(message)

    plasmid_file = str(argv[0])
    common_primer_file = "common_primer.mfasta"
    special_features = "tags_epitopes.mfasta"

    plasmid_records = list(SeqIO.parse(plasmid_file, "genbank"))
    common_primer_records = list(SeqIO.parse(common_primer_file, "fasta"))
    special_features_records = list(SeqIO.parse(special_features, "fasta"))
    
    CommonFeatures().extractFeatures(plasmid_records)
    PrimerBindingSites().extractFeatures(plasmid_records, common_primer_records)

if __name__ == "__main__":
    main(sys.argv[1:])












