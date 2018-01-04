"""
@author: Kevin Streiter & Andreas Ott
"""
from Annotator import *
from Bio import Seq
from Bio.Blast import NCBIWWW
from Bio import SearchIO

class Blaster:
    @staticmethod
    def execute(start, plasmid, blast_direction):
            promotor = False
            end_pos = 0
            addlength = 3 - len(plasmid)%3
            plasmid_to_check = plasmid.seq + plasmid.seq[:3+addlength]

            for i in range(start, len(plasmid_to_check), 3):
                if promotor == False:
                    if plasmid_to_check[i:i + 3] == "ATG":
                        start_pos = i
                        promotor = True
                if promotor == True:
                    if plasmid_to_check[i:i + 3] == "TAG" or plasmid_to_check[i:i + 3] == "TAA" or plasmid_to_check[i:i + 3] == "TGA":
                        OUTPUT_BLAST = "OutputFiles/blasts.xml"
                        end_pos = i
                        end_pos = Annotator().evaluateEndPosition(end_pos, plasmid, len(plasmid.seq))
                        promotor = False
                        translated_hit = Seq.translate(plasmid_to_check[start_pos:end_pos])
                        if len(translated_hit) > 50:
                            print ("Found: " + str(start_pos) + " to " + str(end_pos) + " Size " + str(end_pos - start_pos))
                            blasts = NCBIWWW.qblast("blastp", "nr", str(translated_hit), hitlist_size=1,expect=0.0001)
                            blast_output = open(OUTPUT_BLAST, "w+")
                            blast_output.write(blasts.read())
                            blast_output.close()
                            blast_results = SearchIO.read(OUTPUT_BLAST, 'blast-xml')
                            print (blast_results)
                            print ("")
                            if len(blast_results.hits) > 0:
                                description = blast_results[0].description
                                qualifier = {"note": description, "translation": translated_hit}
                                Annotator().appendFeatures(plasmid, start_pos, end_pos, blast_direction, "CDS", qualifier, 'join')
    @staticmethod
    def blastSearch(plasmid_list):
        for plasmid in plasmid_list:
            Blaster().execute(0,plasmid,1)
            Blaster().execute(0,plasmid.reverse_complement(),-1)
            Blaster().execute(1,plasmid,1)
            Blaster().execute(1,plasmid.reverse_complement(),-1)
            Blaster().execute(2,plasmid,1)
            Blaster().execute(2,plasmid.reverse_complement(),-1)
        Annotator().writeGeneBankFile(plasmid_list, "final_annotated_plasmid.gb")
