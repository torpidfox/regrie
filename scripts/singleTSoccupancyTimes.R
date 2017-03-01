
# collect TS single occupancy times

r45 = read.csv("/Users/dmitrav/Politech/Laboratory/StoсhasticModelling/reGRiE/new/allResults/ts_sim/Raleigh45.csv", sep = " ", header = F)
r850 = read.csv("/Users/dmitrav/Politech/Laboratory/StoсhasticModelling/reGRiE/new/allResults/ts_sim/Raleigh850.csv", sep = " ", header = F)
r336 = read.csv("/Users/dmitrav/Politech/Laboratory/StoсhasticModelling/reGRiE/new/allResults/ts_sim/Raleigh336.csv", sep = " ", header = F)
r391 = read.csv("/Users/dmitrav/Politech/Laboratory/StoсhasticModelling/reGRiE/new/allResults/ts_sim/Raleigh391.csv", sep = " ", header = F)

# find unique TS "ids" in all files
uniqueTS = unique(c(paste0(r45$V1,r45$V2),
                    paste0(r850$V1,r850$V2),
                    paste0(r336$V1,r336$V2),
                    paste0(r391$V1,r391$V2)))

filenames = c()

for (ts in uniqueTS){
    
    time45 = r45[paste0(r45$V1,r45$V2) == ts,"V3"]
    time850 = r850[paste0(r850$V1,r850$V2) == ts,"V3"]
    time336 = r336[paste0(r336$V1,r336$V2) == ts,"V3"]
    time391 = r391[paste0(r391$V1,r391$V2) == ts,"V3"]
    
    maxLength = max(length(time45), length(time850), length(time336), length(time391))
    
    if (maxLength >= 50) {
        
        filenames = c(filenames, paste0(ts,".csv"))
        
        time45 = c(time45, rep(NA, maxLength - length(time45)))
        time850 = c(time850, rep(NA, maxLength - length(time850)))
        time336 = c(time336, rep(NA, maxLength - length(time336)))
        time391 = c(time391, rep(NA, maxLength - length(time391)))
        
        tempdf = data.frame(Raleigh45 = time45,
                            Raleigh850 = time850,
                            Raleigh336 = time336,
                            Raleigh391 = time391)
        
        write.csv(tempdf, paste0("/Users/dmitrav/Politech/Laboratory/StoсhasticModelling/reGRiE/new/allResults/ts_sim/targetSitesByLines/", ts, ".csv"))
    }
    
}

#read files and apply K-S test
pathToFolder = "/Users/dmitrav/Politech/Laboratory/StoсhasticModelling/reGRiE/new/allResults/ts_sim/targetSitesByLines/"

#for testing lines 45 and 850 
p1 = c() # p values
tsid1 = c() # TFname + position

#for testing lines 45 and 850 
p2 = c() # p values
tsid2 = c() # TFname + position

for (fname in filenames){
    
    tempdf = read.csv(paste0(pathToFolder,fname))
    
    test1 = ks.test(tempdf$Raleigh45, tempdf$Raleigh850)
    test2 = ks.test(tempdf$Raleigh336, tempdf$Raleigh391)
    
    if (test1$p.value < 0.001){
        p1 = c(p1, test1$p.value)
        tsid1 = c(tsid1, fname)
    }
    
    if (test2$p.value < 0.001){
        p2 = c(p2, test2$p.value)
        tsid2 = c(tsid2, fname)
    }
    
}

tsid1

