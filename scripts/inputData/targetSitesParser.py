__author__ = 'andrei'

input_file_path = "/Users/dmitrav/Politech/Laboratory/StoсhasticModelling/reGRiE2/biodata/Drosophila/new data/seq_arina_new4.ann"
output_file_path = "/Users/dmitrav/Politech/Laboratory/StoсhasticModelling/reGRiE2/biodata/Drosophila/Kr_TS.txt"

sequence_length = 18001
coding_region_start = 12000
coding_region_end = 14920

gene_name = "Kr"
chromosome_name = "chr2R"
delimiter = ":"

TF_ids = {"hb": 0,"Kr": 1, "gt": 2, "kni": 3, "bcd": 4, "cad": 5, "tll": 6, "hkb": 7}
TF_motif_lengths = {"hb": 10, "Kr": 11, "gt": 12, "kni": 13, "bcd": 7, "cad": 10, "tll": 10, "hkb": 10}


with open(input_file_path, 'r') as file:
    lines = file.readlines()

with open(output_file_path, 'w') as ts:

    for i in range(len(lines)):

        if lines[i].find(">"+gene_name):
            for j in range(1,1+int(lines[i].split(" ")[1])):

                TF_name = lines[j].split(" ")[2].replace("\n","")
                site_end = 1 + sequence_length - int(lines[j].split(" ")[0])
                site_start = site_end - TF_motif_lengths[TF_name]
                DNA_side = "0" if lines[j].split(" ")[1] == '-' else "1"

                if site_end < coding_region_start or site_start > coding_region_start:
                    ts.write(TF_name + delimiter +
                             chromosome_name + delimiter +
                             str(site_start) + ".." +
                             str(site_end) + delimiter +
                             DNA_side + '\n')
            break
