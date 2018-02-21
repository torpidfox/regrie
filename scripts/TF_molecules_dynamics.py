
def plot_TF_dynamics(filename,start_pos_dna,end_pos_dna,first_record_num,num_of_records,plot_id):

    records = []

    # read file and get all data needed
    with open(filename) as file:
        for line in file:

            if line.find('bound at') > 0:
                pos = int(line.split('position')[1].split('in')[0])
                if start_pos_dna <= pos <= end_pos_dna:
                    TF_id = int(line.split('TF')[1].split('of')[0])
                    time = float(line.split(':')[0])
                    TF_pos = pos
                    records.append((TF_id,time,TF_pos))

            if line.find('slide') > 0:
                pos = int(line.split('position')[1].split('by')[0])
                if start_pos_dna <= pos <= end_pos_dna:
                    TF_id = int(line.split('TF')[1].split('of')[0])
                    time = float(line.split(':')[0])
                    TF_pos = pos + 1 if line.find('right') > 0 else pos - 1
                    records.append((TF_id,time,TF_pos))

            if line.find('hop right') > 0 or line.find('hop left') > 0:
                pos = int(line.split('position')[1].split('by')[0])
                if start_pos_dna <= pos <= end_pos_dna:
                    TF_id = int(line.split('TF')[1].split('of')[0])
                    time = float(line.split(':')[0])
                    hop_length = int(line.split('by')[1].split('bp')[0])
                    TF_pos = pos + hop_length if line.find('right') > 0 else pos - hop_length
                    records.append((TF_id, time, TF_pos))

    # sort by time
    import operator
    records = sorted(records,key=operator.itemgetter(1))

    # collecting data for plotting
    labels = []
    datasets = []
    for i in range(first_record_num,num_of_records+first_record_num):
        if records[i][0] in labels:
            datasets[labels.index(records[i][0])][0].append(records[i][1])
            datasets[labels.index(records[i][0])][1].append(records[i][2])
        else:
            labels.append(records[i][0])
            datasets.append(([],[]))
            datasets[labels.index(records[i][0])][0].append(records[i][1])
            datasets[labels.index(records[i][0])][1].append(records[i][2])

    # plotting
    import matplotlib.pyplot as plt

    plt.figure(figsize=(12, 6))
    for i in range(len(labels)):
        plt.plot(datasets[labels.index(labels[i])][0],
                 datasets[labels.index(labels[i])][1],
                 label="TF "+str(labels[i]))

    plt.grid()
    plt.legend(loc=2,prop={'size':10})
    plt.ylabel("Position on the DNA")
    plt.xlabel("Time")
    # plt.show()
    plt.savefig("/Users/dmitrav/Politech/Laboratory/StoсhasticModelling/reGRiE/simulations/status/"+plot_id+".svg")


file = '/Users/dmitrav/Politech/Laboratory/StoсhasticModelling/reGRiE/simulations/status/drosophila_status_3358784414414190651.txt'
plot_TF_dynamics(file,20,120,4500,2000,"8")


