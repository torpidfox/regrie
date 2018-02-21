__author__ = 'andrei'

with open("/home/andrei/Politech/Laboratory/StoсhasticModelling/GRiE-master/biodata/Drosophila/sites_in_output.txt", 'r') as file:
    lines = file.readlines()

hbEnergies = []
KrEnergies = []
gtEnergies = []
kniEnergies = []
bcdEnergies = []
cadEnergies = []
tllEnergies = []
hkbEnergies = []

for line in lines:
    if len(line.split('\t')) > 4:
        if line.split('\t')[3] == '0':
            hbEnergies.append(line.split('\t')[4])
        if line.split('\t')[3] == '1':
            KrEnergies.append(line.split('\t')[4])
        if line.split('\t')[3] == '2':
            gtEnergies.append(line.split('\t')[4])
        if line.split('\t')[3] == '3':
            kniEnergies.append(line.split('\t')[4])
        if line.split('\t')[3] == '4':
            bcdEnergies.append(line.split('\t')[4])
        if line.split('\t')[3] == '5':
            cadEnergies.append(line.split('\t')[4])
        if line.split('\t')[3] == '6':
            tllEnergies.append(line.split('\t')[4])
        if line.split('\t')[3] == '7':
            hkbEnergies.append(line.split('\t')[4])

bcdEnergies = [x for x in bcdEnergies if float(x) > 0]
cadEnergies = [x for x in cadEnergies if float(x) > 0]
tllEnergies = [x for x in tllEnergies if float(x) > 0]

print("hb min energy " + min(hbEnergies))
print("Kr min energy " + min(KrEnergies))
print("gt min energy " + min(gtEnergies))
print("kni min energy " + min(kniEnergies))
print("bcd min energy " + min(bcdEnergies))
print("cad min energy " + min(cadEnergies))
print("tll min energy " + min(tllEnergies))
print("hkb min energy " + min(hkbEnergies))