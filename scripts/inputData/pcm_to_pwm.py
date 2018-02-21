
import os, numpy

pwm = {}

background_probs = {'A': 0.28768819278776,
                    'C': 0.21231180721224,
                    'G': 0.21231180721224,
                    'T': 0.28768819278776}

path = "/Users/dmitrav/Politech/Laboratory/StoсhasticModelling/reGRiE2/debug/pcm/"
for file in os.listdir(path):

    gene = file.split(".")[0]


    pwm[gene] = {'A': [], 'C': [], 'G': [], 'T': []}

    with open(path+file) as f:
        lines = f.readlines()

    for line in lines:

        line = line.replace('\n','').strip()

        if line[0] not in '<>':

            pcm_values = [float(x) for x in line.split(" ")]
            weight = sum(pcm_values)

            pwm[gene]['A'].append(numpy.log((pcm_values[0] + numpy.log(weight) * background_probs['A']) /
                                            (weight + numpy.log(weight)) / background_probs['A']))
            pwm[gene]['C'].append(numpy.log((pcm_values[1] + numpy.log(weight) * background_probs['C']) /
                                            (weight + numpy.log(weight)) / background_probs['C']))
            pwm[gene]['G'].append(numpy.log((pcm_values[2] + numpy.log(weight) * background_probs['G']) /
                                            (weight + numpy.log(weight)) / background_probs['G']))
            pwm[gene]['T'].append(numpy.log((pcm_values[3] + numpy.log(weight) * background_probs['T']) /
                                            (weight + numpy.log(weight)) / background_probs['T']))

for nucleotide in pwm['Kr']:
    print(nucleotide, pwm['Kr'][nucleotide])
    print()
