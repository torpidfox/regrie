
# path to 5e4 simulations series
path = "/Users/dmitrav/Politech/Laboratory/StoсhasticModelling/reGRiE/versions/v.1.2/simulations/5e4/"

folder = c()
for (i in 1:6){
    folder = c(folder, paste("set",i,"-5/results/", sep = ''))
}

# read 5e4 simulations series
occupancy10 = read.csv(paste(path, folder[1], "drosophila_occupancy_1451412716765524957.wig", sep = ''), skip = 1)
occupancy20 = read.csv(paste(path, folder[2], "drosophila_occupancy_5555800597190648070.wig", sep = ''), skip = 1)
occupancy30 = read.csv(paste(path, folder[3], "drosophila_occupancy_3832086354561422017.wig", sep = ''), skip = 1)
occupancy40 = read.csv(paste(path, folder[4], "drosophila_occupancy_146273115075827645.wig", sep = ''), skip = 1)
occupancy50 = read.csv(paste(path, folder[5], "drosophila_occupancy_5343998565444059093.wig", sep = ''), skip = 1)
occupancy60 = read.csv(paste(path, folder[6], "drosophila_occupancy_8159506971092418224.wig", sep = ''), skip = 1)

occupancyList = list(occupancy10, occupancy20, occupancy30, occupancy40, occupancy50, occupancy60)

for (k in 1:length(occupancyList)){
    
    # occupancy 5' -> 3'
    sum53 <-c()
    for (j in 1:nrow(occupancyList[[k]])){
        sum = 0
        for (i in c(3,5,7,9,11,13,15,17)){
            sum <- sum + occupancyList[[k]][j,i]
        }
        sum53 <- c(sum53, sum)
    }

    # occupancy 3' -> 5'
    sum35 <-c()
    for (j in 1:nrow(occupancyList[[k]])){
        sum = 0
        for (i in c(2,4,6,8,10,12,14,16,18)){
            sum <- sum + occupancyList[[k]][j,i]
        }
        sum35 <- c(sum35, sum)
    }

    #add values
    occupancyList[[k]]$seconds_5_3 <- sum53 / 10**6 #in seconds
    occupancyList[[k]]$seconds_3_5 <- sum35 / 10**6 #in seconds
}

#density plots

# 10 seconds
ggplot(data = occupancyList[[1]], aes(occupancyList[[1]]$seconds_5_3)) +
    geom_histogram(aes(y =..density..),alpha = .15,binwidth = 0.001) + 
    geom_density(col=1) +
    labs(title="Density plot for DNA occupancy (10 sec, 5e4 simulations)") +
    labs(x="Time occupied, sec", y="Frequency") + 
    scale_x_continuous(limits = c(0, 0.45)) + 
    scale_y_continuous(breaks = seq(0,1400,100))

# 20 seconds
ggplot(data = occupancyList[[2]], aes(occupancyList[[2]]$seconds_5_3)) +
    geom_histogram(aes(y =..density..),alpha = .15,binwidth = 0.001) + 
    geom_density(col=1) +
    labs(title="Density plot for DNA occupancy (20 sec, 5e4 simulations)") +
    labs(x="Time occupied, sec", y="Frequency") + 
    scale_x_continuous(limits = c(0, 0.45)) + 
    scale_y_continuous(limits = c(0,500), breaks = seq(0,500,100))

# 30 seconds
ggplot(data = occupancyList[[3]], aes(occupancyList[[3]]$seconds_5_3)) +
    geom_histogram(aes(y =..density..),alpha = .15,binwidth = 0.001) + 
    geom_density(col=1) +
    labs(title="Density plot for DNA occupancy (30 sec, 5e4 simulations)") +
    labs(x="Time occupied, sec", y="Frequency") + 
    scale_x_continuous(limits = c(0, 0.45)) + 
    scale_y_continuous(limits = c(0,500), breaks = seq(0,500,100))

# 40 seconds
ggplot(data = occupancyList[[4]], aes(occupancyList[[4]]$seconds_5_3)) +
    geom_histogram(aes(y =..density..),alpha = .15,binwidth = 0.001) + 
    geom_density(col=1) +
    labs(title="Density plot for DNA occupancy (40 sec, 5e4 simulations)") +
    labs(x="Time occupied, sec", y="Frequency") + 
    scale_x_continuous(limits = c(0, 0.45)) + 
    scale_y_continuous(limits = c(0,500), breaks = seq(0,500,100))

# 50 seconds
ggplot(data = occupancyList[[5]], aes(occupancyList[[5]]$seconds_5_3)) +
    geom_histogram(aes(y =..density..),alpha = .15,binwidth = 0.001) + 
    geom_density(col=1) +
    labs(title="Density plot for DNA occupancy (50 sec, 5e4 simulations)") +
    labs(x="Time occupied, sec", y="Frequency") + 
    scale_x_continuous(limits = c(0, 0.45)) + 
    scale_y_continuous(limits = c(0,500), breaks = seq(0,500,100))

# 60 seconds
ggplot(data = occupancyList[[6]], aes(occupancyList[[6]]$seconds_5_3)) +
    geom_histogram(aes(y =..density..),alpha = .15,binwidth = 0.001) + 
    geom_density(col=1) +
    labs(title="Density plot for DNA occupancy (60 sec, 5e4 simulations)") +
    labs(x="Time occupied, sec", y="Frequency") + 
    scale_x_continuous(limits = c(0, 0.45)) + 
    scale_y_continuous(limits = c(0,500), breaks = seq(0,500,100))

#___________________________________________________________________________________________
#___________________________________________________________________________________________

# ABOUT TO RUN ALL THE LINES BELOW B-)

# path to 5e4 simulations series
path2 = "/Users/dmitrav/Politech/Laboratory/StoсhasticModelling/reGRiE/versions/v.1.2/simulations/5e3/"

folder2 = c()
for (i in 1:10){
    folder2 = c(folder2, paste("set",i,"/results/", sep = ''))
}

occupancy1 = read.csv(paste(path2, folder2[1], "drosophila_occupancy_1903667753342154032.wig", sep = ''), skip = 1)
occupancy2 = read.csv(paste(path2, folder2[2], "drosophila_occupancy_1107091672295493714.wig", sep = ''), skip = 1)
occupancy3 = read.csv(paste(path2, folder2[3], "drosophila_occupancy_2089276637961256992.wig", sep = ''), skip = 1)
occupancy4 = read.csv(paste(path2, folder2[4], "drosophila_occupancy_8358990977366430726.wig", sep = ''), skip = 1)
occupancy5 = read.csv(paste(path2, folder2[5], "drosophila_occupancy_7815981773192780744.wig", sep = ''), skip = 1)
occupancy6 = read.csv(paste(path2, folder2[6], "drosophila_occupancy_4137898804458427859.wig", sep = ''), skip = 1)
occupancy7 = read.csv(paste(path2, folder2[7], "drosophila_occupancy_4978610342156924963.wig", sep = ''), skip = 1)
occupancy8 = read.csv(paste(path2, folder2[8], "drosophila_occupancy_8784500193092677272.wig", sep = ''), skip = 1)
occupancy9 = read.csv(paste(path2, folder2[9], "drosophila_occupancy_8462998575499707150.wig", sep = ''), skip = 1)
occupancy10_2 = read.csv(paste(path2, folder2[10], "drosophila_occupancy_3276797477763947272.wig", sep = ''), skip = 1)

occupancyList2 = list(occupancy1, occupancy2, occupancy3, occupancy4, occupancy5, occupancy6, occupancy7, occupancy8, occupancy9, occupancy10_2)

for (k in 1:length(occupancyList2)){
    
    # occupancy 5' -> 3'
    sum53 <-c()
    for (j in 1:nrow(occupancyList2[[k]])){
        sum = 0
        for (i in c(3,5,7,9,11,13,15,17)){
            sum <- sum + occupancyList2[[k]][j,i]
        }
        sum53 <- c(sum53, sum)
    }
    
    # occupancy 3' -> 5'
    sum35 <-c()
    for (j in 1:nrow(occupancyList2[[k]])){
        sum = 0
        for (i in c(2,4,6,8,10,12,14,16,18)){
            sum <- sum + occupancyList2[[k]][j,i]
        }
        sum35 <- c(sum35, sum)
    }
    
    #add values
    occupancyList2[[k]]$seconds_5_3 <- sum53 / 10**6 #in seconds
    occupancyList2[[k]]$seconds_3_5 <- sum35 / 10**6 #in seconds
}

#density plots

# 1 second
ggplot(data = occupancyList2[[1]], aes(occupancyList2[[1]]$seconds_5_3)) +
    geom_histogram(aes(y =..density..),alpha = .15,binwidth = 0.001) + 
    geom_density(col=1) +
    labs(title="Density plot for DNA occupancy (1 sec, 5e3 simulations)") +
    labs(x="Time occupied, sec", y="Frequency") + 
    scale_x_continuous(limits = c(0, 0.1))

# 2 seconds
ggplot(data = occupancyList2[[2]], aes(occupancyList2[[2]]$seconds_5_3)) +
    geom_histogram(aes(y =..density..),alpha = .15,binwidth = 0.001) + 
    geom_density(col=1) +
    labs(title="Density plot for DNA occupancy (2 sec, 5e3 simulations)") +
    labs(x="Time occupied, sec", y="Frequency") + 
    scale_x_continuous(limits = c(0, 0.175))

# 3 seconds
ggplot(data = occupancyList2[[3]], aes(occupancyList2[[3]]$seconds_5_3)) +
    geom_histogram(aes(y =..density..),alpha = .15,binwidth = 0.001) + 
    geom_density(col=1) +
    labs(title="Density plot for DNA occupancy (3 sec, 5e3 simulations)") +
    labs(x="Time occupied, sec", y="Frequency") + 
    scale_x_continuous(limits = c(0, 0.25))

# 4 seconds
ggplot(data = occupancyList2[[4]], aes(occupancyList2[[4]]$seconds_5_3)) +
    geom_histogram(aes(y =..density..),alpha = .15,binwidth = 0.001) + 
    geom_density(col=1) +
    labs(title="Density plot for DNA occupancy (4 sec, 5e3 simulations)") +
    labs(x="Time occupied, sec", y="Frequency") + 
    scale_x_continuous(limits = c(0, 0.375))

# 5 seconds
ggplot(data = occupancyList2[[5]], aes(occupancyList2[[5]]$seconds_5_3)) +
    geom_histogram(aes(y =..density..),alpha = .15,binwidth = 0.001) + 
    geom_density(col=1) +
    labs(title="Density plot for DNA occupancy (5 sec, 5e3 simulations)") +
    labs(x="Time occupied, sec", y="Frequency") + 
    scale_x_continuous(limits = c(0, 0.45))

# 6 seconds
ggplot(data = occupancyList2[[6]], aes(occupancyList2[[6]]$seconds_5_3)) +
    geom_histogram(aes(y =..density..),alpha = .15,binwidth = 0.001) + 
    geom_density(col=1) +
    labs(title="Density plot for DNA occupancy (6 sec, 5e3 simulations)") +
    labs(x="Time occupied, sec", y="Frequency") + 
    scale_x_continuous(limits = c(0, 0.45))

# 7 seconds
ggplot(data = occupancyList2[[7]], aes(occupancyList2[[7]]$seconds_5_3)) +
    geom_histogram(aes(y =..density..),alpha = .15,binwidth = 0.001) + 
    geom_density(col=1) +
    labs(title="Density plot for DNA occupancy (7 sec, 5e3 simulations)") +
    labs(x="Time occupied, sec", y="Frequency") + 
    scale_x_continuous(limits = c(0, 0.45))

# 8 seconds
ggplot(data = occupancyList2[[8]], aes(occupancyList2[[8]]$seconds_5_3)) +
    geom_histogram(aes(y =..density..),alpha = .15,binwidth = 0.001) + 
    geom_density(col=1) +
    labs(title="Density plot for DNA occupancy (8 sec, 5e3 simulations)") +
    labs(x="Time occupied, sec", y="Frequency") + 
    scale_x_continuous(limits = c(0, 0.45))

# 9 seconds
ggplot(data = occupancyList2[[9]], aes(occupancyList2[[9]]$seconds_5_3)) +
    geom_histogram(aes(y =..density..),alpha = .15,binwidth = 0.001) + 
    geom_density(col=1) +
    labs(title="Density plot for DNA occupancy (9 sec, 5e3 simulations)") +
    labs(x="Time occupied, sec", y="Frequency") + 
    scale_x_continuous(limits = c(0, 0.45))

# 10-2 seconds
ggplot(data = occupancyList2[[10]], aes(occupancyList2[[10]]$seconds_5_3)) +
    geom_histogram(aes(y =..density..),alpha = .15,binwidth = 0.001) + 
    geom_density(col=1) +
    labs(title="Density plot for DNA occupancy (10 sec, 5e3 simulations)") +
    labs(x="Time occupied, sec", y="Frequency") + 
    scale_x_continuous(limits = c(0, 0.45))