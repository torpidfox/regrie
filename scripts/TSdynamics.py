import csv, operator
import matplotlib.pyplot as plt

def getStuff(filename):
    with open(filename, "rt", encoding="utf-8") as csvfile:
        datareader = csv.reader(csvfile)
        for row in datareader:
            yield row


# method reads huge file line by line, computes ratio and appends to a list
# after list is large enough it's written to another file and cleared
def processHugeFileForRatio(path,inFileName,outFileName):

    print("Start.\n")
    points = [(0.0, 0.0)]  # [ (time, ratio) ]
    for row in getStuff(path+inFileName):
        ratio = 0
        if row[0].strip() != "time":
            for i in range(1,len(row)):
                ratio += int(row[i])
            if float(row[0].strip()) > 0:
                points.append((float(row[0].strip()),float(ratio)/(len(row)-1)))
        if len(points) == 10 ** 6:
            print("Printing to file...\n")
            with open(outFileName, 'a') as file:
                for point in points:
                    file.write(str(point[0]) + " " + str(point[1]) + "\n")
            points.clear()
            print("List cleared. Collecting data...\n")

    if len(points) > 0:
        with open(outFileName, 'a') as file:
            for point in points:
                file.write(str(point[0]) + " " + str(point[1]) + "\n")
            print("End.")

# method read NOT LARGE file, computes ratio,
# appends to a list of tuple, sorts it by time and plots
def getRatioPlot(path, inFileName):

    print("Start.\n")
    points = [(0.0, 0.0)]  # [ (time, ratio) ]
    for row in getStuff(path + inFileName):
        ratio = 0
        if row[0].strip() != "time":
            for i in range(1, len(row)):
                ratio += int(row[i])
            if float(row[0].strip()) > 0:
                points.append((float(row[0].strip()), float(ratio) / (len(row) - 1)))

    print("Processing...")
    x_points = []
    y_points = []
    for point in sorted(points,key=operator.itemgetter(0)):
        x_points.append(point[0])
        y_points.append(point[1])

    print("Plotting...")
    plt.plot(x_points,y_points,label="Ratio", lw=0.5)
    plt.legend()
    plt.savefig(path+"ratio.svg")
    print("Plot saved.")

# method reads file line by step, get processed values, appends to a list,
# sorts the list by time and plots
def getRatioPlotFromProcessedFile(path, fileName, step):

    points = []
    count = 0
    print("Collecting data...")
    for row in getStuff(path+fileName):
        if count == step:
            points.append((float(row[0].strip().split(" ")[0]), float(row[0].strip().split(" ")[1])))
            count = 0
        count += 1

    x_points = []
    y_points = []
    print("Processing...")
    for point in sorted(points, key=operator.itemgetter(0)):
        x_points.append(point[0])
        y_points.append(point[1])

    print("Plotting...")
    plt.plot(x_points, y_points, label="Ratio", lw=0.5, color="black")
    plt.legend()
    plt.title("Отношение числа занятых сайтов к их общему числу")
    plt.xlabel("Время, с")
    plt.savefig(path + "ratio.svg")
    print("Plot saved.")

def getSingleSiteDynamicsFile(path,inFileName,outFileName,TSofInterestLogic,intermediaryPrintStep):

    if isinstance(TSofInterestLogic,list):
        print("Start.\n")
        indexOfTSofInterest = [-1 for x in TSofInterestLogic]
        TSrecords = [[] for x in TSofInterestLogic]
        time = []
        for row in getStuff(path + inFileName):
            if row[0].strip() == "time":
                for i in range(1,len(row)):
                    if row[i].replace("\"","").strip() in TSofInterestLogic:
                        indexOfTSofInterest[TSofInterestLogic.index(row[i].replace("\"","").strip())] = i
            else:
                for i in range(len(indexOfTSofInterest)):
                    TSrecords[i].append(row[indexOfTSofInterest[i]].strip())
                time.append(row[0].strip())

            if len(time) * len(TSrecords) >= intermediaryPrintStep:
                print("Printing to file...\n")
                with open(path + outFileName, 'a') as file:
                    for i in range(len(TSrecords[0])):
                        for j in range(len(TSrecords)):
                            file.write(TSrecords[j][i] + ", ")
                        file.write(time[i] + "\n")

                # clear all
                time.clear()
                for records in TSrecords:
                    records.clear()
                print("Lists cleared. Collecting data...\n")

        # final records
        if len(time) > 0:
                with open(path + outFileName, 'a') as file:
                    for i in range(len(TSrecords[0])):
                        for j in range(len(TSrecords)):
                            file.write(TSrecords[j][i] + ", ")
                        file.write(time[i] + "\n")
                print("End.")

def getPlotOfSingleTSdynamics(path, fileName, numberOfPointsToPlot, numberOfSites, TSlabels):

    records = [[] for i in range(numberOfSites)]
    # getting data row by step
    print("Collecting data...")
    for row in getStuff(path+fileName):
        for i in range(numberOfSites):
            records[i].append((float(row[numberOfSites]),float(row[i])))

    # sorting by time
    print("Processing...")
    for i in range(numberOfSites):
        records[i] = sorted(records[i],key=operator.itemgetter(0))

    # plotting
    print("Plotting...")
    fig, axes = plt.subplots(numberOfSites,1,sharex=False)

    t = [[] for i in range(numberOfSites)]
    y = [[] for i in range(numberOfSites)]

    import random
    index = random.randint(0, len(records[0]))
    for i in range(numberOfSites):

        if index+numberOfPointsToPlot < len(records[i]):
            for j in range(index,index+numberOfPointsToPlot):
                t[i].append(records[i][j][0])
                y[i].append(records[i][j][1])
        else:
            index -= numberOfPointsToPlot
            for j in range(index, index+numberOfPointsToPlot):
                t[i].append(records[i][j][0])
                y[i].append(records[i][j][1])

        if i == 0: color = 'red'
        elif i == 1: color = 'green'
        else: color = 'blue'

        axes[i].fill_between(t[i], 0, y[i], facecolors=color, alpha=1)
        axes[i].set_ylabel(TSlabels[i],labelpad=20)
    plt.setp(axes, yticks=[0, 1])

    print("Lengths:")
    for i in range(len(t)):
        print(len(t[i]), len(y[i]))

    import numpy
    print("Identity:")
    if len(y[0]) == len(y[1]) == len(y[2]):
        print(numpy.sum(y[0] == y[1] == y[2]) / len(y[0]))
    else:
        print("None")

    plt.tight_layout()
    plt.savefig(path + "singleTSdynamics"+str(numberOfPointsToPlot)+"points1.svg")
    # plt.show()
    print("Plot saved.")


def getSingleTSratioPlot(path,fileName,numberOfSites,TSlabels):

    records = [[] for i in range(numberOfSites)]
    # getting data row by step
    print("Collecting data...")
    for row in getStuff(path + fileName):
        for i in range(numberOfSites):
            records[i].append((float(row[numberOfSites]), float(row[i])))

    # sorting by time
    print("Processing...")
    for i in range(numberOfSites):
        records[i] = sorted(records[i], key=operator.itemgetter(0))

    t = []
    ratio = [[] for i in range(numberOfSites)]

    print("Plotting...")
    for i in range(numberOfSites):
        timeOccupied = 0.
        for j in range(1,len(records[i])):
            if int(records[i][j-1][1]) == 1:
                timeOccupied += records[i][j][0]-records[i][j-1][0]
            ratio[i].append(timeOccupied/records[i][j][0])
            t.append(records[i][j][0])

        # plotting
        plt.figure()
        plt.plot(t, ratio[i], label=TSlabels[i], lw=0.5)
        plt.title("Single site occupancy time ratio")
        plt.legend()
        plt.savefig(path+"singleTSratio"+str(i+1)+".svg" )

        ratio[i].clear()
        t.clear()

    print("All plots saved.")


if __name__ == '__main__':

    # path = "/Volumes/Seagate Expansion Drive/simulations/reGRiE2/Kr_60sec_repression_off/"
    # path2 = "/Users/dmitrav/Politech/Laboratory/StoсhasticModelling/reGRiE/simulations/"
    # inputFile = "drosophila_target_site_follow_884799113612498116.csv"
    # outputFile = "singleSiteDynamics30sec.txt"

    # TS = ["hb:chr3L:14239..14249:0","Kr:chr3L:6140..6151:1","cad:chr3L:15490..15500:1"]
    # getSingleSiteDynamicsFile(path,inputFile,outputFile, TS, 10**6)
    # getPlotOfSingleTSdynamics(path2,outputFile,100,len(TS),[ts.replace(":chr3L:",":") for ts in TS])

    # TS = ["hb:chr3L:14239..14249:0","hb:chr3L:8421..8431:0","hb:chr3L:13596..13606:1"]
    # getSingleSiteDynamicsFile(path,inputFile,outputFile, TS, 10**6)
    # getSingleTSratioPlot(path,outputFile,len(TS),[ts.replace(":chr3L:",":") for ts in TS])

    # TS = ["testSite1","testSite2"]
    # path = "/Users/dmitrav/Downloads/Telegram Desktop/"
    # outputFile = "testDynamics.txt"
    # getSingleTSratioPlot(path,outputFile,len(TS),[ts for ts in TS])

    path = "/Users/dmitrav/Politech/Laboratory/StoсhasticModelling/reGRiE2/simulations/Kr_60sec_repression_off/"
    # path2 = "/Users/dmitrav/Politech/Laboratory/StoсhasticModelling/reGRiE/simulations/"
    inputFile = "drosophila_target_site_follow_884799113612498116.csv"
    outputFile = "ratio_data.txt"

    # processHugeFileForRatio(path,inputFile,path+outputFile)
    getRatioPlotFromProcessedFile(path,outputFile,5)





