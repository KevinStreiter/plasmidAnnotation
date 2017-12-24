"""
@author: Kevin Streiter & Andreas Ott
"""
import sys
from CommonFeatures import *
from Bio import SeqIO


def main(argv):

    if len(argv) == 0:
        message = "Please provide a file"
        sys.exit(message)

    filename = str(argv[0])
    plasmid_records = SeqIO.parse(filename, "genbank")

    CommonFeatures(plasmid_records)

if __name__ == "__main__":
    main(sys.argv[1:])












