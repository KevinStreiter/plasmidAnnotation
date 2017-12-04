# -*- coding: utf-8 -*-
"""
@author: Kevin Streiter
"""
from Bio import SeqIO

def scanFile(filename):
    numbers = []
    with open(filename) as f:
        numbers.append(zip(*[line.split() for line in f])[1])
    return list(map(int, numbers[0]))

def saveGenes(record):
    genes = []
    for gene in record.features:
        if(gene.type == "CDS"):
            genes.append(gene)
    return genes

def index_genbank_features(genes, qualifier):
    answer = dict()
    for (index, feature) in enumerate(genes):
        if qualifier in feature.qualifiers:
            for value in feature.qualifiers[qualifier]:
                if value in answer :
                    print "WARNING - Duplicate key, please specify query"
                    return answer
                else:
                    answer[value] = index
    return answer


def calculateNumberOfSequencedBases(start, end, numbers):
    return sum(numbers[start:end])

cellulose = 'Cellulose__GSM463010_080606_HWI-EAS282_0007_FC305F9AAXX.7_coverage.txt'
glukose = 'Glukose__GSM463006_080527_HWI-EAS282_0004_FC30316AAXX.1_coverage.txt'

cellulose_numbers = scanFile(cellulose)
glukose_numbers = scanFile(glukose)

print "--------------NUMBER OF SEQUENCED BASES--------------"
print "Cellulose:", sum( cellulose_numbers)
print "Glukose:", sum(glukose_numbers)

record = SeqIO.read(open("genome_sso.gb","r"),"genbank")
genes = saveGenes(record)
sequence = record.seq
translated_record = sequence.reverse_complement().translate(table=11)

locus_tag_index = index_genbank_features(genes, "locus_tag")
index = locus_tag_index["SSO0564"]

gene = genes[index]

start = gene.location.start.position-1
end = gene.location.end.position-1

print calculateNumberOfSequencedBases(start, end, cellulose_numbers)
print calculateNumberOfSequencedBases(start, end, glukose_numbers)










