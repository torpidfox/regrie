with open("/home/andrei/Politech/Laboratory/StoсhasticModelling/GRiE-master/biodata/Drosophila/sites_in_output.txt", 'r') as file:
    lines = file.readlines()

hbExpEnergies = []
KrExpEnergies = []
gtExpEnergies = []
kniExpEnergies = []
bcdExpEnergies = []
cadExpEnergies = []
tllExpEnergies = []
hkbExpEnergies = []

for line in lines:
    if len(line.split('\t')) > 4:
        if line.split('\t')[3] == '0':
            hbExpEnergies.append(line.split('\t')[5])
        if line.split('\t')[3] == '1':
            KrExpEnergies.append(line.split('\t')[5])
        if line.split('\t')[3] == '2':
            gtExpEnergies.append(line.split('\t')[5])
        if line.split('\t')[3] == '3':
            kniExpEnergies.append(line.split('\t')[5])
        if line.split('\t')[3] == '4':
            bcdExpEnergies.append(line.split('\t')[5])
        if line.split('\t')[3] == '5':
            cadExpEnergies.append(line.split('\t')[5])
        if line.split('\t')[3] == '6':
            tllExpEnergies.append(line.split('\t')[5])
        if line.split('\t')[3] == '7':
            hkbExpEnergies.append(line.split('\t')[5])

hbAverageExp = 0
for value in hbExpEnergies:
    hbAverageExp += float(value)
hbAverageExp /= len(hbExpEnergies)

KrAverageExp = 0
for value in KrExpEnergies:
    KrAverageExp += float(value)
KrAverageExp /= len(KrExpEnergies)

gtAverageExp = 0
for value in gtExpEnergies:
    gtAverageExp += float(value)
gtAverageExp /= len(gtExpEnergies)

kniAverageExp = 0
for value in kniExpEnergies:
    kniAverageExp += float(value)
kniAverageExp /= len(kniExpEnergies)

bcdAverageExp = 0
for value in bcdExpEnergies:
    bcdAverageExp += float(value)
bcdAverageExp /= len(bcdExpEnergies)

cadAverageExp = 0
for value in cadExpEnergies:
    cadAverageExp += float(value)
cadAverageExp /= len(cadExpEnergies)

tllAverageExp = 0
for value in tllExpEnergies:
    tllAverageExp += float(value)
tllAverageExp /= len(tllExpEnergies)

hkbAverageExp = 0
for value in hkbExpEnergies:
    hkbAverageExp += float(value)
hkbAverageExp /= len(hkbExpEnergies)

print("hb average exp(energy) = " + str(hbAverageExp))
print("Kr average exp(energy) = " + str(KrAverageExp))
print("gt average exp(energy) = " + str(gtAverageExp))
print("kni average exp(energy) = " + str(kniAverageExp))
print("bcd average exp(energy) = " + str(bcdAverageExp))
print("cad average exp(energy) = " + str(cadAverageExp))
print("tll average exp(energy) = " + str(tllAverageExp))
print("hkb average exp(energy) = " + str(hkbAverageExp))

#calculate specific waiting time according to Zabet 2012
#tau_0 = 2 * t_R / s_1_obs ^ 2 / mean(exp(PWM))
with open("specificWaitingTime.txt", "w") as file:
    file.write("hb tau_0 = " + str(round(4 / 425 / 425 / hbAverageExp, 7)) + "\n")
    file.write("Kr tau_0 = " + str(round(4 / 425 ** 2 / KrAverageExp, 7)) + "\n")
    file.write("gt tau_0 = " + str(round(4 / 425 ** 2 / gtAverageExp, 7)) + "\n")
    file.write("kni tau_0 = " + str(round(4 / 425 ** 2 / kniAverageExp, 7)) + "\n")
    file.write("bcd tau_0 = " + str(round(4 / 425 ** 2 / bcdAverageExp, 7)) + "\n")
    file.write("cad tau_0 = " + str(round(4 / 425 ** 2 / cadAverageExp, 7)) + "\n")
    file.write("tll tau_0 = " + str(round(4 / 425 ** 2 / tllAverageExp, 7)) + "\n")
    file.write("hkb tau_0 = " + str(round(4 / 425 ** 2 / hkbAverageExp, 7)) + "\n")