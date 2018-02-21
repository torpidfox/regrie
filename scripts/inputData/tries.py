with open("/Users/dmitrav/Politech/Laboratory/StoсhasticModelling/reGRiE2/biodata/Drosophila/genes&chromatin/gt_Raleigh_45.fasta") as file:
    seq = file.readlines()[1]

print(seq.find("ACCGAAAATC"))
print(seq[3246:3253])
print(seq[11958:11958+10])

print(seq[7646:7656])