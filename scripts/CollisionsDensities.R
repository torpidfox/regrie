
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

#density plots

# 10 seconds
    #collisions
ggplot(data = occupancyList[[1]], aes(occupancyList[[1]]$collisionsCount)) +
    geom_histogram(aes(y =..density..),alpha = .15,binwidth = 12) + 
    geom_density(col=1) +
    labs(title="Density plot for collisions (10 sec, 5e4 simulations)") +
    labs(x="Number of collisions", y="Frequency") + 
    scale_x_continuous(limits = c(0, 1000)) + 
    scale_y_continuous(limits = c(0.000,0.016),breaks = seq(0.000,0.016,0.002))

# 20 seconds
    #collisions
ggplot(data = occupancyList[[2]], aes(occupancyList[[2]]$collisionsCount)) +
    geom_histogram(aes(y =..density..),alpha = .15,binwidth = 12) + 
    geom_density(col=1) +
    labs(title="Density plot for collisions (20 sec, 5e4 simulations)") +
    labs(x="Number of collisions", y="Frequency") + 
    scale_x_continuous(limits = c(0, 1000)) + 
    scale_y_continuous(limits = c(0.000,0.016),breaks = seq(0.000,0.016,0.002))

# 30 seconds
    #collisions
ggplot(data = occupancyList[[3]], aes(occupancyList[[3]]$collisionsCount)) +
    geom_histogram(aes(y =..density..),alpha = .15,binwidth = 12) + 
    geom_density(col=1) +
    labs(title="Density plot for collisions (30 sec, 5e4 simulations)") +
    labs(x="Number of collisions", y="Frequency") + 
    scale_x_continuous(limits = c(0, 1000)) +
    scale_y_continuous(limits = c(0.000,0.016),breaks = seq(0.000,0.016,0.002))

# 40 seconds
    #collisions
ggplot(data = occupancyList[[4]], aes(occupancyList[[4]]$collisionsCount)) +
    geom_histogram(aes(y =..density..),alpha = .15,binwidth = 12) + 
    geom_density(col=1) +
    labs(title="Density plot for collisions (40 sec, 5e4 simulations)") +
    labs(x="Number of collisions", y="Frequency") + 
    scale_x_continuous(limits = c(0, 1000)) +
    scale_y_continuous(limits = c(0.000,0.016),breaks = seq(0.000,0.016,0.002))

# 50 seconds
    #collisions
ggplot(data = occupancyList[[5]], aes(occupancyList[[5]]$collisionsCount)) +
    geom_histogram(aes(y =..density..),alpha = .15,binwidth = 12) + 
    geom_density(col=1) +
    labs(title="Density plot for collisions (50 sec, 5e4 simulations)") +
    labs(x="Number of collisions", y="Frequency") + 
    scale_x_continuous(limits = c(0, 1000)) +
    scale_y_continuous(limits = c(0.000,0.016),breaks = seq(0.000,0.016,0.002))

# 60 seconds
    #collisions
ggplot(data = occupancyList[[6]], aes(occupancyList[[6]]$collisionsCount)) +
    geom_histogram(aes(y =..density..),alpha = .15,binwidth = 12) + 
    geom_density(col=1) +
    labs(title="Density plot for collisions (60 sec, 5e4 simulations)") +
    labs(x="Number of collisions", y="Frequency") + 
    scale_x_continuous(limits = c(0, 1000)) +
    scale_y_continuous(limits = c(0.000,0.016),breaks = seq(0.000,0.016,0.002))

df = data.frame(collisions = occupancy10$collisionsCount, seconds = 10)
df30 = data.frame(collisions = occupancy30$collisionsCount, seconds = 30)
df60 = data.frame(collisions = occupancy60$collisionsCount, seconds = 60)

total = rbind(df, df30, df60)

ggplot(data = total, aes(x = collisions, fill = as.factor(seconds))) +
    geom_density(alpha = .15) +
    labs(x="Число столкновений", y="Частота") + 
    scale_x_continuous(limits = c(0, 1000)) +
    scale_y_continuous(limits = c(0.000,0.016),breaks = seq(0.000,0.016,0.002)) +
    scale_fill_discrete(name="Время симуляций, с") + 
    theme(plot.title = element_text(vjust=20))


set.seed(1234)
dat <- data.frame(cond = factor(rep(c("A","B"), each=200)), 
                  rating = c(rnorm(200),rnorm(200, mean=.8)))

ggplot(dat, aes(x=rating, fill=cond)) + geom_density(alpha=.3)
