__author__ = 'andrei'

with open("/Users/dmitrav/Politech/Laboratory/StoсhasticModelling/reGRiE/versions/reGRiE v.1.2/biodata/Drosophila/sites_in_output.txt", 'r') as file:
    lines = file.readlines()

with open("/Users/dmitrav/Politech/Laboratory/StoсhasticModelling/reGRiE/versions/reGRiE v.1.2/biodata/Drosophila/ts.csv", 'w') as ts:

    for line in lines:
        if len(line.split('\t')) > 4:

            if line.split('\t')[3] == '0' \
                    and (int(line.split('\t')[1]) + int(line.split('\t')[7]) < 4144 or int(line.split('\t')[1]) > 6001) \
                    and line.split('\t')[2] == '-':
                ts.write("hb:chr3L:" + line.split('\t')[1] + ".." + str((int(line.split('\t')[1]) + int(line.split('\t')[7]))) + ":0" + '\n')

            if line.split('\t')[3] == '0' \
                    and (int(line.split('\t')[1]) + int(line.split('\t')[7]) < 4144 or int(line.split('\t')[1]) > 6001) \
                    and line.split('\t')[2] == '+':
                ts.write("hb:chr3L:" + line.split('\t')[1] + ".." + str((int(line.split('\t')[1]) + int(line.split('\t')[7]))) + ":1" + '\n')

            if line.split('\t')[3] == '1' \
                    and (int(line.split('\t')[1]) + int(line.split('\t')[7]) < 4144 or int(line.split('\t')[1]) > 6001) \
                    and line.split('\t')[2] == '-':
                ts.write("Kr:chr3L:" + line.split('\t')[1] + ".." + str((int(line.split('\t')[1]) + int(line.split('\t')[7]))) + ":0" + '\n')

            if line.split('\t')[3] == '1' \
                    and (int(line.split('\t')[1]) + int(line.split('\t')[7]) < 4144 or int(line.split('\t')[1]) > 6001) \
                    and line.split('\t')[2] == '+':
                ts.write("Kr:chr3L:" + line.split('\t')[1] + ".." + str((int(line.split('\t')[1]) + int(line.split('\t')[7]))) + ":1" + '\n')

            if line.split('\t')[3] == '2' \
                    and (int(line.split('\t')[1]) + int(line.split('\t')[7]) < 4144 or int(line.split('\t')[1]) > 6001) \
                    and line.split('\t')[2] == '-':
                ts.write("gt:chr3L:" + line.split('\t')[1] + ".." + str((int(line.split('\t')[1]) + int(line.split('\t')[7]))) + ":0" + '\n')

            if line.split('\t')[3] == '2' \
                    and (int(line.split('\t')[1]) + int(line.split('\t')[7]) < 4144 or int(line.split('\t')[1]) > 6001) \
                    and line.split('\t')[2] == '+':
                ts.write("gt:chr3L:" + line.split('\t')[1] + ".." + str((int(line.split('\t')[1]) + int(line.split('\t')[7]))) + ":1" + '\n')

            if line.split('\t')[3] == '3' \
                    and (int(line.split('\t')[1]) + int(line.split('\t')[7]) < 4144 or int(line.split('\t')[1]) > 6001) \
                    and line.split('\t')[2] == '-':
                ts.write("kni:chr3L:" + line.split('\t')[1] + ".." + str((int(line.split('\t')[1]) + int(line.split('\t')[7]))) + ":0" + '\n')

            if line.split('\t')[3] == '3' \
                    and (int(line.split('\t')[1]) + int(line.split('\t')[7]) < 4144 or int(line.split('\t')[1]) > 6001) \
                    and line.split('\t')[2] == '+':
                ts.write("kni:chr3L:" + line.split('\t')[1] + ".." + str((int(line.split('\t')[1]) + int(line.split('\t')[7]))) + ":1" + '\n')

            if line.split('\t')[3] == '4' \
                    and (int(line.split('\t')[1]) + int(line.split('\t')[7]) < 4144 or int(line.split('\t')[1]) > 6001) \
                    and line.split('\t')[2] == '-':
                ts.write("bcd:chr3L:" + line.split('\t')[1] + ".." + str((int(line.split('\t')[1]) + int(line.split('\t')[7]))) + ":0" + '\n')

            if line.split('\t')[3] == '4' \
                    and (int(line.split('\t')[1]) + int(line.split('\t')[7]) < 4144 or int(line.split('\t')[1]) > 6001) \
                    and line.split('\t')[2] == '+':
                ts.write("bcd:chr3L:" + line.split('\t')[1] + ".." + str((int(line.split('\t')[1]) + int(line.split('\t')[7]))) + ":1" + '\n')

            if line.split('\t')[3] == '5' \
                    and (int(line.split('\t')[1]) + int(line.split('\t')[7]) < 4144 or int(line.split('\t')[1]) > 6001) \
                    and line.split('\t')[2] == '-':
                ts.write("cad:chr3L:" + line.split('\t')[1] + ".." + str((int(line.split('\t')[1]) + int(line.split('\t')[7]))) + ":0" + '\n')

            if line.split('\t')[3] == '5' \
                    and (int(line.split('\t')[1]) + int(line.split('\t')[7]) < 4144 or int(line.split('\t')[1]) > 6001) \
                    and line.split('\t')[2] == '+':
                ts.write("cad:chr3L:" + line.split('\t')[1] + ".." + str((int(line.split('\t')[1]) + int(line.split('\t')[7]))) + ":1" + '\n')

            if line.split('\t')[3] == '6' \
                    and (int(line.split('\t')[1]) + int(line.split('\t')[7]) < 4144 or int(line.split('\t')[1]) > 6001) \
                    and line.split('\t')[2] == '-':
                ts.write("tll:chr3L:" + line.split('\t')[1] + ".." + str((int(line.split('\t')[1]) + int(line.split('\t')[7]))) + ":0" + '\n')

            if line.split('\t')[3] == '6' \
                    and (int(line.split('\t')[1]) + int(line.split('\t')[7]) < 4144 or int(line.split('\t')[1]) > 6001) \
                    and line.split('\t')[2] == '+':
                ts.write("tll:chr3L:" + line.split('\t')[1] + ".." + str((int(line.split('\t')[1]) + int(line.split('\t')[7]))) + ":1" + '\n')

            if line.split('\t')[3] == '7' \
                    and (int(line.split('\t')[1]) + int(line.split('\t')[7]) < 4144 or int(line.split('\t')[1]) > 6001) \
                    and line.split('\t')[2] == '-':
                ts.write("hkb:chr3L:" + line.split('\t')[1] + ".." + str((int(line.split('\t')[1]) + int(line.split('\t')[7]))) + ":0" + '\n')

            if line.split('\t')[3] == '7' \
                    and (int(line.split('\t')[1]) + int(line.split('\t')[7]) < 4144 or int(line.split('\t')[1]) > 6001) \
                    and line.split('\t')[2] == '+':
                ts.write("hkb:chr3L:" + line.split('\t')[1] + ".." + str((int(line.split('\t')[1]) + int(line.split('\t')[7]))) + ":1" + '\n')