
occupancyFile = "/Users/dmitrav/Politech/Laboratory/StoсhasticModelling/reGRiE/versions/v.1.2/simulations/5e4/set1-5/results/drosophila_occupancy_1451412716765524957.wig"
TSfile = "/Users/dmitrav/Politech/Laboratory/StoсhasticModelling/reGRiE/versions/v.1.2/biodata/Drosophila/TF_Drosophila_ts2.csv"
TSinfoFile = "/Users/dmitrav/Politech/Laboratory/StoсhasticModelling/reGRiE/versions/v.1.2/biodata/Drosophila/copy_sites_in_output.txt"

#read files
occupancy = read.csv(occupancyFile, skip = 1)
ts_original = read.csv(TSfile, sep = " ", header = FALSE, skip = 1)
ts = read.csv(TSfile, sep = " ", header = FALSE, skip = 1)
tsIO = read.csv(TSinfoFile, sep = "\t")

#convert indices to names
for (i in 1:(length(tsIO$Factor.Idx))){
    if (tsIO$Factor.Idx[i] == 0){ tsIO$Factor.Idx[i] = 'hb' }
    if (tsIO$Factor.Idx[i] == 1){ tsIO$Factor.Idx[i] = 'Kr' }
    if (tsIO$Factor.Idx[i] == 2){ tsIO$Factor.Idx[i] = 'gt' }
    if (tsIO$Factor.Idx[i] == 3){ tsIO$Factor.Idx[i] = 'kni' }
    if (tsIO$Factor.Idx[i] == 4){ tsIO$Factor.Idx[i] = 'bcd' }
    if (tsIO$Factor.Idx[i] == 5){ tsIO$Factor.Idx[i] = 'cad' }
    if (tsIO$Factor.Idx[i] == 6){ tsIO$Factor.Idx[i] = 'tll' }
    if (tsIO$Factor.Idx[i] == 7){ tsIO$Factor.Idx[i] = 'hkb' }
}

#take energies
ts$E = 0
for (i in 1:length(ts$V1)){
    for (j in 1:length(tsIO$Factor.Idx)){
        if (ts$V1[i] == tsIO$Factor.Idx[j] & ts$V3[i] == tsIO$Start...1[j]){ 
            ts$E[i] = tsIO$Energy[j]
            break
        }
    }
}

#compute exponents
ts$expE = exp(-ts$E)

#available positions for giant enhancer
availabilityBoundaries <- c(3024, 3279, 6002, 6319, 6965, 8713, 8890, 9196, 9914, 10321, 12015, 12335, 13541, 14291, 15265, 16290, 16712, 17079)
              
#recompute TS positions
for (j in 1:nrow(ts)){
    sum <- 0
    i <- 1
    while (ts[j,"V4"] > availabilityBoundaries[i+1]){
        sum <- sum + availabilityBoundaries[i+1] - availabilityBoundaries[i];
        i <- i + 2;
    }
    ts[j,"V3"] <- sum + ts[j,"V3"] - availabilityBoundaries[i];
    ts[j,"V4"] <- sum + ts[j,"V4"] - availabilityBoundaries[i];
}

#create data frame
tsdf = data.frame(name = ts$V1, position = ts$V3, direction = ts$V5, pwm = ts$E, energy = ts$expE, time = 0)

#fill with values for each TF in each orientation

TFnames = c("tll", "Kr", "hb", "cad", "hkb", "gt", "kni", "bcd")

for (name in TFnames){
    for (i in 1:length(tsdf$position)){
        if (tsdf$name[i] == name){
            if (tsdf$direction[i] == 1){
                tsdf$time[i] = occupancy[occupancy$position == tsdf$position[i],paste(name,"5.3.", sep = "")]
            }
            if (tsdf$direction[i] == 0){
                tsdf$time[i] = occupancy[occupancy$position == tsdf$position[i],paste(name,"3.5.", sep = "")]
            }
        }
    }
}

#plotting by pwm
ggplot(tsdf, aes(x = pwm)) +
    geom_line(aes(y = time / 10 ** 6)) +
    geom_hline(aes(yintercept = 30000 / 10 ** 6), colour = 2) + 
    scale_y_continuous(limits = c(0,0.035), breaks = seq(0,0.035,0.005)) + 
    ylab("Time occupied, seconds") +
    xlab("TS Energies") +
    labs(title = "10 seconds simulation")

#plotting by exp(pwm)
ggplot(tsdf, aes(x = energy)) +
    geom_line(aes(y = time / 10 ** 6)) +
    geom_hline(aes(yintercept = 30000 / 10 ** 6), colour = 2) + 
    scale_y_continuous(limits = c(0,0.035), breaks = seq(0,0.035,0.005)) + 
    ylab("Time occupied, seconds") +
    xlab("TS Energies") +
    labs(title = "10 seconds simulation")


