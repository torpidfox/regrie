with open("/Users/dmitrav/Politech/Laboratory/StoсhasticModelling/reGRiE2/biodata/Drosophila/new data/seqNmod4.fa") \
    as file:
    sequences_data = file.readlines()

gene_of_interest = "Kr"

Kr_sequence = sequences_data[1+sequences_data.index(">"+gene_of_interest+"\n")].strip()

with open("/Users/dmitrav/Politech/Laboratory/StoсhasticModelling/reGRiE2/biodata/Drosophila/genes&chromatin/Kr_dnaseAccS05.btrack") \
    as another_file:
    data = another_file.readlines()

btrack_data = ""

for value in data:
    btrack_data += value.replace("\n","")

print(Kr_sequence)
print(btrack_data)

new_btrack_data = ""
for i in range(len(Kr_sequence)):
    if Kr_sequence[i] == "N":
        new_btrack_data += "0"
    else:
        new_btrack_data += "1"

print(new_btrack_data)

# with open("/Users/dmitrav/Politech/Laboratory/StoсhasticModelling/reGRiE2/biodata/Drosophila/genes&chromatin/Kr_corrected.btrack",'w') \
#     as new_btrack_file:
#
#     for value in new_btrack_data:
#         new_btrack_file.write(value+"\n")
