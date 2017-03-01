
simTSpath = "/Users/dmitrav/Politech/Laboratory/StoсhasticModelling/reGRiE/new/allResults/Raleigh850/results/Raleigh_850/drosophila_target_site_8678851071760239436.csv"
gtTSpath = "/Users/dmitrav/Politech/Laboratory/StoсhasticModelling/reGRiE/new/allResults/Raleigh45/biodata/Drosophila/gtTS_space.csv"
SIOpath = "/Users/dmitrav/Politech/Laboratory/StoсhasticModelling/reGRiE/new/allResults/Raleigh45/biodata/Drosophila/new/sites_in_output.txt"

simTS = read.csv(simTSpath)
gtTSspace = read.csv(gtTSpath, sep = " ", header = F)
sio = read.csv(SIOpath, sep = '\t')

#rename sio TFs
for (i in 1:(length(sio$Factor.Idx))){
    if (sio$Factor.Idx[i] == 0){ sio$Factor.Idx[i] = 'hb' }
    if (sio$Factor.Idx[i] == 1){ sio$Factor.Idx[i] = 'Kr' }
    if (sio$Factor.Idx[i] == 2){ sio$Factor.Idx[i] = 'gt' }
    if (sio$Factor.Idx[i] == 3){ sio$Factor.Idx[i] = 'kni' }
    if (sio$Factor.Idx[i] == 4){ sio$Factor.Idx[i] = 'bcd' }
    if (sio$Factor.Idx[i] == 5){ sio$Factor.Idx[i] = 'cad' }
    if (sio$Factor.Idx[i] == 6){ sio$Factor.Idx[i] = 'tll' }
    if (sio$Factor.Idx[i] == 7){ sio$Factor.Idx[i] = 'hkb' }
}

#copy energies to gtTSspace
for (tfName in unique(sio$Factor.Idx)){
    tfdf = subset(sio, sio$Factor.Idx == tfName)
    for (i in 1:nrow(tfdf)){
        gtTSspace[gtTSspace$V1 == tfName & gtTSspace$V3 == (18002 - tfdf$Start...1[i] - tfdf$length[i]),"Energy"] = tfdf$Energy[i]
    }
}

#concatenate TFs' info to get a string of 'logic'
gtTSspace$logic = paste0(gtTSspace$V1,":",gtTSspace$V2,":",gtTSspace$V3,"..",gtTSspace$V4,":",gtTSspace$V5)

#filter simTS
#simTS = simTS[simTS$firstReached > 0,]
simTS[simTS$timeOccupied < 0, "timeOccupied"] = abs(simTS[simTS$timeOccupied < 0, "timeOccupied"])

output = cbind(output,Raleigh850 = simTS$timeOccupied)
output[,c(2,3,4,5)] = output[,c(2,3,4,5)] / 10

output = data.frame(TS = paste0(gtTSspace$V1,":chrX:",gtTSspace$V3,"..",gtTSspace$V4,":",gtTSspace$V5), Raleigh391 = simTS$timeOccupied)


write.csv(output, "/Users/dmitrav/Politech/Laboratory/StoсhasticModelling/reGRiE/new/allResults/TSaverageTime.csv")

#add energies to simTS
for (i in 1:nrow(simTS)){
    simTS$Energy[i] = gtTSspace[simTS$targetSite[i] == gtTSspace$logic,"Energy"]
}

#plotting
ggplot(simTS, aes(x = Energy)) + 
    geom_line(aes(y = timeOccupied / 10)) +
    ylab("Время занятости, мс") + 
    xlab("Энергия сайта") + 
    scale_x_continuous(breaks = seq(0,8,1)) + 
    scale_y_continuous(breaks = seq(0,350,50)) +
    ggtitle("Время занятости целевых сайтов линии Raleigh 391")

ggplot(data = simTS, aes(simTS$timeOccupied / 10)) +
    geom_histogram(aes(y =..density..),alpha = .1,binwidth = 10) + 
    geom_density(col=1) +
    labs(x="Время занятости, мс", y="Частота") +
    scale_x_continuous(limits = c(0,400), breaks = seq(0,400,50)) +
    scale_y_continuous(limits = c(0,0.08)) +
    ggtitle("Частота времен занятости целевых сайтов линии Raleigh 391")

#plotting
ggplot(output, aes(x = Energy)) + 
    geom_line(aes(y = Raleigh850)) +
    ylab("Time occupied, ms") + 
    xlab("Site energy, PWM score") + 
    scale_x_continuous(breaks = seq(0,8,1)) + 
    scale_y_continuous(breaks = seq(0,350,50)) +
    ggtitle("Raleigh 850")

ggplot(data = output, aes(output$Raleigh850)) +
    geom_histogram(aes(y =..density..),alpha = .1,binwidth = 10) + 
    geom_density(col=1) +
    labs(x="Time occupied, ms", y="Frequency") +
    scale_x_continuous(limits = c(0,400), breaks = seq(0,400,50)) +
    scale_y_continuous(limits = c(0,0.08)) +
    ggtitle("Raleigh 850")


