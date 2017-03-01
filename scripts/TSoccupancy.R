path = "/Users/dmitrav/Politech/Laboratory/StoсhasticModelling/reGRiE/versions/v.1.2/simulations/5e3/"
path2 = "/Users/dmitrav/Politech/Laboratory/StoсhasticModelling/reGRiE/versions/v.1.2/simulations/5e4/"

folder = c()
for (i in 1:10){
    folder = c(folder, paste("set",i,"/results/", sep = ''))
}

occupancy1 = read.csv(paste(path, folder[1], "drosophila_occupancy_1903667753342154032.wig", sep = ''), skip = 1)
occupancy2 = read.csv(paste(path, folder[2], "drosophila_occupancy_1107091672295493714.wig", sep = ''), skip = 1)
occupancy3 = read.csv(paste(path, folder[3], "drosophila_occupancy_2089276637961256992.wig", sep = ''), skip = 1)
occupancy4 = read.csv(paste(path, folder[4], "drosophila_occupancy_8358990977366430726.wig", sep = ''), skip = 1)
occupancy5 = read.csv(paste(path, folder[5], "drosophila_occupancy_7815981773192780744.wig", sep = ''), skip = 1)
occupancy6 = read.csv(paste(path, folder[6], "drosophila_occupancy_4137898804458427859.wig", sep = ''), skip = 1)
occupancy7 = read.csv(paste(path, folder[7], "drosophila_occupancy_4978610342156924963.wig", sep = ''), skip = 1)
occupancy8 = read.csv(paste(path, folder[8], "drosophila_occupancy_8784500193092677272.wig", sep = ''), skip = 1)
occupancy9 = read.csv(paste(path, folder[9], "drosophila_occupancy_8462998575499707150.wig", sep = ''), skip = 1)
occupancy10 = read.csv(paste(path, folder[10], "drosophila_occupancy_3276797477763947272.wig", sep = ''), skip = 1)

occupancy60 = read.csv(paste(path2, "set6-5/results/", "drosophila_occupancy_8159506971092418224.wig", sep = ''), skip = 1)


#read file with same TS separated by spaces
ts <- read.csv("/Users/dmitrav/Politech/Laboratory/StoсhasticModelling/reGRiE/versions/v.1.2/biodata/Drosophila/TF_Drosophila_ts2.csv", sep = " ", header = FALSE)

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

#separate TFs
tll <- ts[ts["V1"] == "tll",]
Kr <- ts[ts["V1"] == "Kr",]
bcd <- ts[ts["V1"] == "bcd",]
kni <- ts[ts["V1"] == "kni",]
hb <- ts[ts["V1"] == "hb",]
gt <- ts[ts["V1"] == "gt",]
cad <- ts[ts["V1"] == "cad",]
hkb <- ts[ts["V1"] == "hkb",]

#hb TS occupancy analysis
#considering length of hb TS of value 10 

#10 seconds
occupancy_sum_10 <- 0
for (i in 1:nrow(hb)){
    for (j in -5:5){
        occupancy_sum_10 <- occupancy_sum_10 + occupancy10[hb[i,"V3"] + j,"hb5.3."] + occupancy10[hb[i,"V3"] + j,"hb3.5."]
    }
}

#20 seconds
occupancy_sum_20 <- 0
for (i in 1:nrow(hb)){
    for (j in -5:5){
        occupancy_sum_20 <- occupancy_sum_20 + occupancy20[hb[i,"V3"] + j,"hb5.3."] + occupancy20[hb[i,"V3"] + j,"hb3.5."]
    }
}

#30 seconds
occupancy_sum_30 <- 0
for (i in 1:nrow(hb)){
    for (j in -5:5){
        occupancy_sum_30 <- occupancy_sum_30 + occupancy30[hb[i,"V3"] + j,"hb5.3."] + occupancy30[hb[i,"V3"] + j,"hb3.5."]
    }
}

#40 seconds
occupancy_sum_40 <- 0
for (i in 1:nrow(hb)){
    for (j in -5:5){
        occupancy_sum_40 <- occupancy_sum_40 + occupancy40[hb[i,"V3"] + j,"hb5.3."] + occupancy40[hb[i,"V3"] + j,"hb3.5."]
    }
}

#50 seconds
occupancy_sum_50 <- 0
for (i in 1:nrow(hb)){
    for (j in -5:5){
        occupancy_sum_50 <- occupancy_sum_50 + occupancy50[hb[i,"V3"] + j,"hb5.3."] + occupancy50[hb[i,"V3"] + j,"hb3.5."]
    }
}

#60 seconds
occupancy_sum_60 <- 0
for (i in 1:nrow(hb)){
    for (j in -5:5){
        occupancy_sum_60 <- occupancy_sum_60 + occupancy60[hb[i,"V3"] + j,"hb5.3."] + occupancy60[hb[i,"V3"] + j,"hb3.5."]
    }
}

hb_TS_occupancy <- data.frame(c(10,20,30,40,50,60), c(occupancy_sum_10, occupancy_sum_20, occupancy_sum_30, occupancy_sum_40, occupancy_sum_50, occupancy_sum_60))
hb_TS_occupancy["scaled_occupancy"] <- scale(hb_TS_occupancy$c.occupancy_sum_10..occupancy_sum_20..occupancy_sum_30..occupancy_sum_40..)

#normalize
hb_TS_occupancy["normalized_occupancy"] <- 0

for (i in 1:nrow(hb_TS_occupancy)){
    hb_TS_occupancy[i,"normalized_occupancy"] <- hb_TS_occupancy[i,"c.occupancy_sum_10..occupancy_sum_20..occupancy_sum_30..occupancy_sum_40.."] / (i * 10**7)
}

ggplot(hb_TS_occupancy, aes(x=c.10..20..30..40..50..60., y=normalized_occupancy)) +
    geom_bar(stat="identity", fill = "grey", colour = "black", width=8) +
    ylab("Frequency of occupied target sites") +
    xlab('Biological time, sec') +
    scale_x_continuous(breaks = seq(0,60,10))

#TS occupancy for all TS
TS_ordered <- ts[order(ts[,3]),]

#Calculate TS occupancy for different time periods
TS_occupancy1 <- c()
TS_occupancy2 <- c()
TS_occupancy3 <- c()
TS_occupancy4 <- c()
TS_occupancy5 <- c()
TS_occupancy6 <- c()
TS_occupancy7 <- c()
TS_occupancy8 <- c()
TS_occupancy9 <- c()
TS_occupancy10 <- c()

TS_occupancy60 <- c()

for (i in 1:nrow(TS_ordered)){
    sum1 <- 0
    sum2 <- 0
    sum3 <- 0
    sum4 <- 0
    sum5 <- 0
    sum6 <- 0
    sum7 <- 0
    sum8 <- 0
    sum9 <- 0
    sum10 <- 0
    
    sum60 = 0
    
    for (j in TS_ordered[i,"V3"]:TS_ordered[i,"V4"]){
        sum1 <- sum1 + occupancy1[j-1,paste(TS_ordered[i,1], "5.3.", sep='')] + occupancy1[j-1,paste(TS_ordered[i,1], "3.5.", sep='')]
        sum2 <- sum2 + occupancy2[j-1,paste(TS_ordered[i,1], "5.3.", sep='')] + occupancy2[j-1,paste(TS_ordered[i,1], "3.5.", sep='')]
        sum3 <- sum3 + occupancy3[j-1,paste(TS_ordered[i,1], "5.3.", sep='')] + occupancy3[j-1,paste(TS_ordered[i,1], "3.5.", sep='')]
        sum4 <- sum4 + occupancy4[j-1,paste(TS_ordered[i,1], "5.3.", sep='')] + occupancy4[j-1,paste(TS_ordered[i,1], "3.5.", sep='')]
        sum5 <- sum5 + occupancy5[j-1,paste(TS_ordered[i,1], "5.3.", sep='')] + occupancy5[j-1,paste(TS_ordered[i,1], "3.5.", sep='')]
        sum6 <- sum6 + occupancy6[j-1,paste(TS_ordered[i,1], "5.3.", sep='')] + occupancy6[j-1,paste(TS_ordered[i,1], "3.5.", sep='')]
        sum7 <- sum7 + occupancy7[j-1,paste(TS_ordered[i,1], "5.3.", sep='')] + occupancy7[j-1,paste(TS_ordered[i,1], "3.5.", sep='')]
        sum8 <- sum8 + occupancy8[j-1,paste(TS_ordered[i,1], "5.3.", sep='')] + occupancy8[j-1,paste(TS_ordered[i,1], "3.5.", sep='')]
        sum9 <- sum9 + occupancy9[j-1,paste(TS_ordered[i,1], "5.3.", sep='')] + occupancy9[j-1,paste(TS_ordered[i,1], "3.5.", sep='')]
        sum10 <- sum10 + occupancy10[j-1,paste(TS_ordered[i,1], "5.3.", sep='')] + occupancy10[j-1,paste(TS_ordered[i,1], "3.5.", sep='')]
    
        sum60 <- sum60 + occupancy60[j-1,paste(TS_ordered[i,1], "5.3.", sep='')] + occupancy60[j-1,paste(TS_ordered[i,1], "3.5.", sep='')]
        
        }
    
    TS_occupancy1 <- c(TS_occupancy1, sum1)
    TS_occupancy2 <- c(TS_occupancy2, sum2)
    TS_occupancy3 <- c(TS_occupancy3, sum3)
    TS_occupancy4 <- c(TS_occupancy4, sum4)
    TS_occupancy5 <- c(TS_occupancy5, sum5)
    TS_occupancy6 <- c(TS_occupancy6, sum6)
    TS_occupancy7 <- c(TS_occupancy7, sum7)
    TS_occupancy8 <- c(TS_occupancy8, sum8)
    TS_occupancy9 <- c(TS_occupancy9, sum9)
    TS_occupancy10 <- c(TS_occupancy10, sum10)
    
    TS_occupancy60 <- c(TS_occupancy60, sum60)
}

TS_ordered$occupancy1 <- TS_occupancy1
TS_ordered$occupancy2 <- TS_occupancy2
TS_ordered$occupancy3 <- TS_occupancy3
TS_ordered$occupancy4 <- TS_occupancy4
TS_ordered$occupancy5 <- TS_occupancy5
TS_ordered$occupancy6 <- TS_occupancy6
TS_ordered$occupancy7 <- TS_occupancy7
TS_ordered$occupancy8 <- TS_occupancy8
TS_ordered$occupancy9 <- TS_occupancy9
TS_ordered$occupancy10 <- TS_occupancy10

TS_ordered$occupancy60 <- TS_occupancy60

#plotting TS occupancy
ggplot(data = TS_ordered, aes(x=1:nrow(TS_ordered))) +
    geom_line(aes(y = occupancy_10 / 10**6, colour = "10 sec")) +
    geom_line(aes(y = occupancy_30 / 10**6, colour = "30 sec")) +
    geom_line(aes(y = occupancy_60 / 10**6, colour ="60 sec")) +
    ylab("Time occupied, sec") +
    xlab("Target sites") +
    scale_x_continuous(breaks = seq(0,260, 20)) +
    scale_colour_manual(values = c("10 sec"="#CF3F42", "30 sec"="#A67800", "60 sec"="#63ADD0"))

#total density plot    
ggplot(data = TS_ordered, aes(TS_ordered$occupancy60 / 10**6)) +
    geom_histogram(aes(y =..density..),alpha = .1,binwidth = 0.08) + 
    geom_density(col=1) +
    labs(x="Время занятости, с", y="Частота")


#density plots    
ggplot(data = TS_ordered, aes(TS_ordered$occupancy10 / 10**6)) +
    geom_histogram(aes(y =..density..),alpha = .3,binwidth = 0.005) + 
    geom_density(col=1) +
    labs(title="Density plot for TS occupancy (10 sec, 5e3 simulations)") +
    labs(x="Time occupied, sec", y="Frequency") + 
    scale_x_continuous(limits = c(0, 0.05))

ggplot(data = TS_ordered, aes(TS_ordered$occupancy_20 / 10**6)) +
    geom_histogram(aes(y =..density..),alpha = .3,binwidth = 0.05) + 
    geom_density(col=1) +
    labs(title="Density plot for TS occupancy (20 seconds SSA)") +
    labs(x="Time occupied, sec", y="Frequency")

ggplot(data = TS_ordered, aes(TS_ordered$occupancy_30 / 10**6)) +
    geom_histogram(aes(y =..density..),alpha = .3,binwidth = 0.05) + 
    geom_density(col=1) +
    labs(title="Density plot for TS occupancy (30 seconds SSA)") +
    labs(x="Time occupied, sec", y="Frequency")

ggplot(data = TS_ordered, aes(TS_ordered$occupancy_40 / 10**6)) +
    geom_histogram(aes(y =..density..),alpha = .3,binwidth = 0.05) + 
    geom_density(col=1) +
    labs(title="Density plot for TS occupancy (40 seconds SSA)") +
    labs(x="Time occupied, sec", y="Frequency")

ggplot(data = TS_ordered, aes(TS_ordered$occupancy_50 / 10**6)) +
    geom_histogram(aes(y =..density..),alpha = .3,binwidth = 0.05) + 
    geom_density(col=1) +
    labs(title="Density plot for TS occupancy (50 seconds SSA)") +
    labs(x="Time occupied, sec", y="Frequency")

ggplot(data = TS_ordered, aes(TS_ordered$occupancy_60 / 10**6)) +
    geom_histogram(aes(y =..density..),alpha = .3,binwidth = 0.05) + 
    geom_density(col=1) +
    labs(title="Density plot for TS occupancy (60 seconds SSA)") +
    labs(x="Time occupied, sec", y="Frequency")


#plays
chol <- read.table(url("http://assets.datacamp.com/blog_assets/chol.txt"), header = TRUE)
ggplot(data=chol, aes(chol$AGE)) + 
    geom_histogram(aes(y =..density..), 
                   breaks=seq(20, 50, by = 2), 
                   col="red", 
                   fill="green", 
                   alpha = .2) + 
    geom_density(col=2) + 
    labs(title="Histogram for Age") +
    labs(x="Age", y="Count")


