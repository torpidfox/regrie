
path = "/Users/dmitrav/Politech/Laboratory/StoсhasticModelling/reGRiE/versions/v.1.2/simulations/5e3/"

folder = "set10/results/"

occupancy = read.csv(paste(path, folder, "drosophila_occupancy_3276797477763947272.wig", sep = ''), skip = 1)
affinity = read.csv(paste(path, folder, "drosophila_affinity_landscape_1903667753342154032.wig", sep = ''), skip = 1)

library(ggplot2)

#collisions

ggplot(data = occupancy, aes(x = position)) +
    geom_bar(aes(y = collisionsCount), stat = "identity") +
    ylab("Collisions count") +
    xlab('Positions on the DNA') +
    scale_x_continuous(breaks = seq(0,nrow(occupancy),300)) + 
    labs(title = "DNA collisions (10 sec, 5e3 simulations)")
    
#occupancy 5' -> 3'
sum53 <-c()
for (j in 1:nrow(occupancy)){
    sum = 0
    for (i in c(3,5,7,9,11,13,15,17)){
        sum <- sum + occupancy[j,i]
    }
    sum53 <- c(sum53, sum)
}

#occupancy 3' -> 5'
sum35 <-c()
for (j in 1:nrow(occupancy)){
    sum = 0
    for (i in c(2,4,6,8,10,12,14,16,18)){
        sum <- sum + occupancy[j,i]
    }
    sum35 <- c(sum35, sum)
}

#add values
occupancy$direction_5_3 <- sum53
occupancy$direction_3_5 <- sum35
occupancy$normalized_5_3 <- sum53 / 10**6 #in seconds
occupancy$normalized_3_5 <- sum35 / 10**6 #in seconds

library(ggplot2)

#summarized occupancy
ggplot(data = occupancy, aes(x = position)) +
    geom_bar(aes(y = normalized_3_5, color = "3' -> 5'"), stat = "identity") +
    geom_bar(aes(y = normalized_5_3, color = "5' -> 3'"), stat = "identity") +
    ylab("Time occupied, sec") +
    xlab('Positions on the DNA') +
    scale_x_continuous(breaks = seq(0,nrow(occupancy),300)) +
    labs(title = "DNA occupancy sum (10 sec, 5e3 simulations)")

#form Kr data
Kr_data <- data.frame(occupancy$Kr5.3., occupancy$Kr3.5., abs(affinity$Kr5.3.), abs(affinity$Kr3.5.))

#Kr occupancy and affinity
ggplot(data = Kr_data, aes(x = 1:nrow(Kr_data))) +
    geom_line(aes(y = occupancy.Kr5.3.)) +
    geom_line(aes(y = abs.affinity.Kr5.3.., color = 'red'))

#scale Kr data
Kr_scaled_occupancy53 <- scale(Kr_data$occupancy.Kr5.3.)
Kr_scaled_affinity53 <- scale(Kr_data$abs.affinity.Kr5.3..)

scaled_Kr_data <- data.frame(Kr_scaled_occupancy53, Kr_scaled_affinity53)

#Kr scaled occupancy and affinity
ggplot(data = scaled_Kr_data, aes(x = 1:nrow(scaled_Kr_data))) +
    geom_line(aes(y = Kr_scaled_occupancy53)) +
    geom_line(aes(y = Kr_scaled_affinity53, color = 'red'))

#form cad data
cad_data <- data.frame(occupancy$cad5.3., occupancy$cad3.5., abs(affinity$cad5.3.), abs(affinity$cad3.5.))

#cad occupancy and affinity
ggplot(data = cad_data, aes(x = 1:nrow(cad_data))) +
    geom_line(aes(y = occupancy.cad5.3.)) +
    geom_line(aes(y = abs.affinity.cad5.3.., color = 'red'))

#scale cad data
cad_scaled_occupancy53 <- scale(cad_data$occupancy.cad5.3.)
cad_scaled_affinity53 <- scale(cad_data$abs.affinity.cad5.3..)

scaled_cad_data <- data.frame(cad_scaled_occupancy53, cad_scaled_affinity53)

#cad scaled occupancy and affinity
ggplot(data = scaled_cad_data, aes(x = 1:nrow(scaled_cad_data))) +
    geom_line(aes(y = cad_scaled_occupancy53)) +
    geom_curve(aes(y = cad_scaled_occupancy53)) +
    geom_line(aes(y = cad_scaled_affinity53, color = 'red'))

