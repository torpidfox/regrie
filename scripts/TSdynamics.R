data = read.csv("/Users/dmitrav/Politech/Laboratory/StoсhasticModelling/reGRiE/versions/v.1.2/results/drosophila_target_site_follow_8747095416411813331.csv", header = TRUE)

#remove sites never occupied
followedTS = data[,which(!apply(data,2,FUN = function(x){all(x == 0)}))]

#don't remove, why remove?
followedTS = data

#select data of ~zero points
zeroState = followedTS[which(followedTS$time == 0 | followedTS$time == 0.000000e+00),]

#time boundary and uniform time intervals
highTime = min(followedTS$time[followedTS$time > 0]) * 10 ** 6
microTime = seq(0,highTime,length.out = length(zeroState$time))

zeroState$time = microTime

#calculate ratio of TS occupancy during time
zeroState$ratio = apply(zeroState[,2:ncol(zeroState)],1,mean)

#plotting
ggplot(zeroState, aes(x = microTime)) + 
    geom_line(aes(y = ratio))

#remove ~zero time points 
followedTS = followedTS[which(followedTS$time > highTime / 10 ** 6),]

#calculate ratio of TS occupancy during time
followedTS$ratio = apply(followedTS[,2:ncol(followedTS)],1,mean)

#ascending sorting by time
followedTS = followedTS[order(followedTS$time),]

#calculate time between two recordings (TBTR)
followedTS$tbtr = 0
for (i in 2:length(followedTS$time)){
    followedTS$tbtr[i] = followedTS$time[i] - followedTS$time[i-1]
}

#plotting
ggplot(followedTS, aes(x = ratio)) + 
    geom_line(aes(y = tbtr))
